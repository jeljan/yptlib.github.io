import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
import html
import py3Dmol
from pathlib import Path
import plotly.express as px
import plotly.graph_objs as go
from shiny import App, ui, render, reactive
from shinywidgets import output_widget, render_widget
from functools import reduce

# --- 1. Setup Paths & Data Loading ---
app_dir = Path(__file__).parent
data_dir = app_dir / "data"

all_files = list(data_dir.glob("*.csv"))
df = None

if len(all_files) > 0:
    df_list = [pd.read_csv(f) for f in all_files]
    df = reduce(lambda left, right: pd.merge(left, right, on='Info', how='outer'), df_list)

    if 'Info' in df.columns:
        df[['Protein Id', 'Gene Symbol', 'Site Position', 'Sequence', 'Description']] = df['Info'].str.split('++', expand=True)
        df.drop(columns='Info', inplace=True)
        df['Site Position'] = df['Site Position'].fillna('Unknown')
        df['Labels'] = df['Gene Symbol'] + "_Y" + df['Site Position'].astype(str)

    raw_drugs = list(set([i.split(' ')[1] for i in df.columns if 'log' in i]))
    
    # Promiscuity Calculation
    drug_promiscuity = {}
    for d in raw_drugs:
        log_col = f'log2 {d} R'
        if log_col in df.columns:
            total_valid_sites = df[log_col].notna().sum()
            if total_valid_sites > 0:
                hits = (df[log_col] > 1).sum()
                promiscuity = (hits / total_valid_sites) * 100
            else:
                promiscuity = 0.0
            drug_promiscuity[d] = promiscuity
            
    sorted_drugs = sorted(raw_drugs, key=lambda x: drug_promiscuity[x], reverse=True)
    drug_choices = {d: f"{d} ({drug_promiscuity[d]:.2f}%)" for d in sorted_drugs}
    default_drug = sorted_drugs[0] if sorted_drugs else None
    
    label_choices = sorted(df['Labels'].dropna().unique().tolist())
else:
    df = pd.DataFrame()
    drug_choices = {"No Data Found": "No Data Found"}
    default_drug = "No Data Found"
    label_choices = ["No Data Found"]

# Safely load the Smiles dictionary & Cancer List
ypt_lib_path = app_dir / 'ypt_library.csv'
if ypt_lib_path.exists():
    ypt_lib = pd.read_csv(ypt_lib_path)
    smiles_dict = dict(zip(ypt_lib['Catalog ID'], ypt_lib['Smiles']))
else:
    smiles_dict = {}
    
cancer_path = app_dir / 'cancer_gene_shortlist.csv'
if cancer_path.exists():
    cancer = pd.read_csv(cancer_path)
else:
    cancer = pd.DataFrame(columns=['Gene']) 

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Tyrosine Library Screening"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_select("data_type", "Select Drug (Promiscuity %):", choices=drug_choices, selected=default_drug),
            ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
            ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
            
            ui.hr(),
            ui.input_select(
                "color_mode", "Color Points By:", 
                choices=["Above Threshold", "P-Value Gradient", "Cancer-Driver List", "Highlight Custom List"]
            ),
            ui.input_text("custom_list", "Genes to Highlight (comma-separated):", placeholder="e.g. MAPK1, EGFR"),
            
            ui.hr(),
            ui.h5("Chemical Structure"),
            ui.output_ui("molecule_ui"), 
            
            width=300
        ),
        
        # --- TABBED LAYOUT ---
        ui.navset_card_tab(
            ui.nav_panel(
                "Compound View", 
                ui.div(
                    ui.h5("Compound Engagement Profile"),
                    output_widget("site_plot"),
                    ui.hr(),
                    ui.h5("Volcano Plot"),
                    output_widget("volcano_plot"),
                    style="padding: 10px; height: 100%; overflow-y: auto;"
                ),
                value="compound_tab"
            ),
            ui.nav_panel(
                "Site View",
                ui.layout_columns(
                    ui.card(
                        ui.h5("Engagement by Compound"),
                        ui.input_selectize("target_site", "Search & Select Site:", choices=label_choices),
                        output_widget("site_profile_plot")
                    ),
                    ui.card(
                        ui.h5("AlphaFold 3D Structure"),
                        ui.p("Selected site is highlighted in red.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                        ui.output_ui("alphafold_viewer")
                    ),
                    col_widths=(5, 7)
                ),
                value="site_tab" 
            ),
            id="main_tabs", 
            full_screen=True
        )
    )
)

# --- 3. Shiny Server ---
def server(input, output, session):
    
    active_drug = reactive.Value(default_drug)
    
    @reactive.Effect
    @reactive.event(input.main_tabs)
    def handle_tab_switch():
        if input.main_tabs() == "site_tab":
            active_drug.set(None) 
        else:
            active_drug.set(input.data_type()) 
            
    @reactive.Effect
    @reactive.event(input.data_type)
    def handle_dropdown():
        if input.main_tabs() == "compound_tab":
            active_drug.set(input.data_type())

    @render.ui
    def molecule_ui():
        drug = active_drug.get()
        if not drug:
            return ui.p("Click a bar in the chart to view structure.", style="color: gray; font-style: italic;")

        smiles = smiles_dict.get(drug, "") 
        if not smiles:
            return ui.p("No SMILES available for this drug.", style="color: gray; font-style: italic;")

        safe_smiles = urllib.parse.quote(smiles)
        img_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/image"
        return ui.HTML(f'<img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Chemical Structure">')

    # --- SHARED REACTIVE DATA FOR BOTH PLOTS ---
    @reactive.Calc
    def plot_data():
        data_type = input.data_type()
        threshold = input.threshold()
        n_labels = input.n_labels()
        color_mode = input.color_mode()
        custom_list_text = input.custom_list()

        plot_df = pd.DataFrame({
            'Site': np.arange(len(df[f'log2 {data_type} R'])),
            'R': 2**df[f'log2 {data_type} R'],
            'p': 10**-df[f'-log10 {data_type} p'],
            'log2 R': df[f'log2 {data_type} R'],
            '-log10 p': df[f'-log10 {data_type} p'],
            'Label': df['Labels'],
            'Gene Symbol': df['Gene Symbol'],
            'ID': df['Protein Id'],
            'Description': df['Description']
        })
        
        plot_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        plot_df.dropna(subset=['R', 'p', 'log2 R', '-log10 p'], inplace=True)

        plot_df['color'] = 'non-significant'

        if color_mode == "Above Threshold":
            plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
            sig_genes = plot_df[plot_df['color'] == 'high']
        elif color_mode == "Highlight Custom List":
            custom_genes = [g.strip().upper() for g in custom_list_text.split(',') if g.strip()]
            mask = plot_df['Gene Symbol'].str.upper().isin(custom_genes)
            plot_df.loc[mask, 'color'] = 'highlight'
            sig_genes = plot_df[plot_df['color'] == 'highlight']
        elif color_mode == "Cancer-Driver List":          
            mask = plot_df['Gene Symbol'].isin(cancer['Gene'])
            plot_df.loc[mask, 'color'] = 'highlight'
            sig_genes = plot_df[plot_df['color'] == 'highlight']
        else: # P-Value Gradient
            plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
            sig_genes = plot_df[plot_df['color'] == 'high']
            plot_df['alpha'] = np.where(plot_df['R'] > threshold, 1.0, 0.2)

        # Label logic
        valid_indices = []
        pos_count = min(n_labels, len(sig_genes))
        if pos_count > 0:
            top_positive = sig_genes.nlargest(pos_count, 'R').index
            valid_indices = [idx for idx in top_positive if idx in plot_df.index]
            
        return plot_df, valid_indices

    # Standard hover formatting dictionary shared between plots
    hover_dict = {
        'Site': False, 'color': False, 
        'ID': True, 'Label': True, 'Description': True,     
        'R': ':.3f', 'p': ':.4f', 'log2 R': ':.3f', '-log10 p': ':.3f'
    }

    # --- TAB 1: PLOT 1 (ENGAGEMENT PLOT) ---
    @render_widget
    def site_plot():
        plot_df, valid_indices = plot_data()
        color_mode = input.color_mode()
        threshold = input.threshold()

        if color_mode == "Above Threshold":
            fig = px.scatter(
                plot_df, x='Site', y='R', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=0.9), selector=dict(name='high'))

        elif color_mode in ["Highlight Custom List", "Cancer-Driver List"]:
            fig = px.scatter(
                plot_df, x='Site', y='R', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'highlight': '#D62728'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.15), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))

        else: # P-Value Gradient
            fig = px.scatter(
                plot_df, x='Site', y='R', color='p',
                color_continuous_scale='Viridis_r', 
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=plot_df['alpha']))

        for idx in valid_indices:
            row = plot_df.loc[idx]
            fig.add_annotation(
                x=row['Site'], y=row['R'], text=row['Label'],
                showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1.5,
                arrowcolor="#444", ax=20, ay=-30,
                font=dict(size=10, color="black"), bgcolor="white",
                bordercolor="#c7c7c7", borderwidth=1, borderpad=3
            )

        fig.add_hline(y=threshold, line_width=1, line_dash='dash')
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))

        return fig

    # --- TAB 1: PLOT 2 (VOLCANO PLOT) ---
    @render_widget
    def volcano_plot():
        plot_df, valid_indices = plot_data()
        color_mode = input.color_mode()
        threshold = input.threshold()

        if color_mode == "Above Threshold":
            fig = px.scatter(
                plot_df, x='log2 R', y='-log10 p', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=0.9), selector=dict(name='high'))

        elif color_mode in ["Highlight Custom List", "Cancer-Driver List"]:
            fig = px.scatter(
                plot_df, x='log2 R', y='-log10 p', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'highlight': '#D62728'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.15), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))

        else: # P-Value Gradient
            fig = px.scatter(
                plot_df, x='log2 R', y='-log10 p', color='p',
                color_continuous_scale='Viridis_r', 
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=plot_df['alpha']))

        for idx in valid_indices:
            row = plot_df.loc[idx]
            fig.add_annotation(
                x=row['log2 R'], y=row['-log10 p'], text=row['Label'],
                showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=1.5,
                arrowcolor="#444", ax=20, ay=-30,
                font=dict(size=10, color="black"), bgcolor="white",
                bordercolor="#c7c7c7", borderwidth=1, borderpad=3
            )

        # Volcano-specific threshold lines
        safe_thresh = threshold if threshold > 0 else 1e-9 # Prevent log2(0) crash
        fig.add_vline(x=np.log2(safe_thresh), line_width=1, line_dash='dash')
        fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash='dash')
        
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))

        return fig


    # --- TAB 2: SITE PROFILING (BAR CHART WITH CLICKS) ---
    @render_widget
    def site_profile_plot():
        target = input.target_site()
        if not target or target == "No Data Found" or df.empty:
            return go.FigureWidget(px.bar(title="No Site Selected"))

        row = df[df['Labels'] == target].iloc[0]
        
        drug_data = []
        for d in raw_drugs: 
            log_col = f'log2 {d} R'
            if log_col in df.columns:
                r_val = 2**row[log_col]
                drug_data.append({'Drug': d, 'R Value': r_val})
                
        bar_df = pd.DataFrame(drug_data).dropna()
        bar_df = bar_df.sort_values('R Value', ascending=False)
        
        fig = px.bar(
            bar_df, x='Drug', y='R Value', 
            title=f"Engagement Profile: {target}", 
            color='R Value', color_continuous_scale='Reds'
        )
        
        fig.add_hline(y=input.threshold(), line_width=1, line_dash='dash', line_color='red', annotation_text="Threshold")
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=20, r=20, t=40, b=20), coloraxis_showscale=False)
        
        widget = go.FigureWidget(fig)
        
        def handle_bar_click(trace, points, state):
            if points.point_inds:
                idx = points.point_inds[0]
                clicked_drug = trace.x[idx]
                active_drug.set(clicked_drug)
                
        if widget.data:
            widget.data[0].on_click(handle_bar_click)
            
        return widget

    # --- TAB 2: ALPHAFOLD VIEWER ---
    @render.ui
    def alphafold_viewer():
        target = input.target_site()
        if not target or target == "No Data Found" or df.empty:
            return ui.p("Please select a site from the dropdown.")

        row = df[df['Labels'] == target].iloc[0]
        
        protein_string = str(row['Protein Id'])
        if '|' in protein_string:
            uniprot_raw = protein_string.split('|')[1]
        else:
            uniprot_raw = protein_string
            
        site_pos = str(row['Site Position'])
        uniprot = uniprot_raw.split('-')[0]

        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.pdb"
        try:
            req = urllib.request.urlopen(url)
            pdb_data = req.read().decode('utf-8')
        except Exception:
            return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve structure for {uniprot}. It may not exist in the AlphaFold database or the network request failed.</div>")

        view = py3Dmol.view(width="100%", height=500)
        view.addModel(pdb_data, "pdb")
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.addStyle({'resi': site_pos}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
        view.zoomTo({'resi': site_pos})
        
        viewer_html = html.escape(view._make_html())
        iframe = f'<iframe srcdoc="{viewer_html}" style="width: 100%; height: 500px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>'
        
        return ui.HTML(iframe)

# --- 4. Run App ---
app = App(app_ui, server)