import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
import html
import py3Dmol
from pathlib import Path
import plotly.express as px
from shiny import App, ui, render
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
        df[['Protein Id', 'Gene Symbol', 'Site Position', 'Sequence']] = df['Info'].str.split('+', expand=True)
        df.drop(columns='Info', inplace=True)
        df['Site Position'] = df['Site Position'].fillna('Unknown')
        df['Labels'] = df['Gene Symbol'] + "_Y" + df['Site Position'].astype(str)

    drugs = list(set([i.split(' ')[1] for i in df.columns if 'log' in i]))
    # Create a clean list of all unique sites for the new dropdown
    label_choices = sorted(df['Labels'].dropna().unique().tolist())
else:
    df = pd.DataFrame()
    drugs = ["No Data Found"]
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
    ui.h2("Proteomics Engagement Dashboard"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.h5("Global Plot Settings"),
            ui.input_select("data_type", "Select Drug:", choices=drugs, selected=drugs[0] if drugs else None),
            ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
            ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
            
            ui.hr(),
            ui.h5("Plot Styling"),
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
        
        # --- NEW TABBED LAYOUT ---
        ui.navset_card_tab(
            # TAB 1: The original volcano-style plot
            ui.nav_panel(
                "Compound Engagement", 
                output_widget("site_plot")
            ),
            
            # TAB 2: The new site-specific visualization and AlphaFold viewer
            ui.nav_panel(
                "Site Profiling & Structure",
                ui.layout_columns(
                    ui.card(
                        ui.h5("Cross-Compound Profiling"),
                        ui.input_selectize("target_site", "Search & Select Site:", choices=label_choices),
                        output_widget("site_profile_plot")
                    ),
                    ui.card(
                        ui.h5("AlphaFold 3D Structure"),
                        ui.p("Selected site is highlighted in red.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                        ui.output_ui("alphafold_viewer")
                    ),
                    col_widths=(5, 7) # Sets the width ratio between the bar chart and 3D viewer
                )
            ),
            full_screen=True
        )
    )
)

# --- 3. Shiny Server ---
def server(input, output, session):
    
    # --- UI RENDERS ---
    @render.ui
    def molecule_ui():
        data_type = input.data_type()
        smiles = smiles_dict.get(data_type, "") 
        if not smiles:
            return ui.p("No SMILES available for this drug.", style="color: gray; font-style: italic;")

        safe_smiles = urllib.parse.quote(smiles)
        img_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/image"
        return ui.HTML(f'<img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Chemical Structure">')

    # --- TAB 1: GLOBAL PLOT ---
    @render_widget
    def site_plot():
        data_type = input.data_type()
        threshold = input.threshold()
        n_labels = input.n_labels()
        color_mode = input.color_mode()
        custom_list_text = input.custom_list()

        plot_df = pd.DataFrame({
            'Site': np.arange(len(df[f'log2 {data_type} R'])),
            'R': 2**df[f'log2 {data_type} R'],
            'p': 10**-df[f'-log10 {data_type} p'],
            'Label': df['Labels'],
            'Gene Symbol': df['Gene Symbol'],
            'ID': df['Protein Id'] + '|' + df['Labels'].str.split('_', expand=True).iloc[:, 1]
        })
        
        plot_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        plot_df.dropna(subset=['R', 'p'], inplace=True)

        plot_df['color'] = 'non-significant'

        if color_mode == "Above Threshold":
            plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
            sig_genes = plot_df[plot_df['color'] == 'high']
            
            fig = px.scatter(
                plot_df, x='Site', y='R', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD'},
                hover_data=['ID', 'R', 'p']
            )
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=0.9), selector=dict(name='high'))

        elif color_mode == "Highlight Custom List":
            custom_genes = [g.strip() for g in custom_list_text.split(',') if g.strip()]
            mask = plot_df['Gene Symbol'].isin(custom_genes)
            plot_df.loc[mask, 'color'] = 'highlight'
            sig_genes = plot_df[plot_df['color'] == 'highlight']
            
            fig = px.scatter(
                plot_df, x='Site', y='R', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'highlight': '#D62728'},
                hover_data=['ID', 'R', 'p']
            )
            fig.update_traces(marker=dict(opacity=0.15), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))
            
        elif color_mode == "Cancer-Driver List":          
            mask = plot_df['Gene Symbol'].isin(cancer['Gene'])
            plot_df.loc[mask, 'color'] = 'highlight'
            sig_genes = plot_df[plot_df['color'] == 'highlight']
            
            fig = px.scatter(
                plot_df, x='Site', y='R', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'highlight': '#D62728'},
                hover_data=['ID', 'R', 'p']
            )
            fig.update_traces(marker=dict(opacity=0.15), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))

        else: # P-Value Gradient
            plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
            sig_genes = plot_df[plot_df['color'] == 'high']
            plot_df['alpha'] = np.where(plot_df['R'] > threshold, 1.0, 0.2)
            
            fig = px.scatter(
                plot_df, x='Site', y='R', color='p',
                color_continuous_scale='Viridis', 
                hover_data=['ID', 'R', 'p']
            )
            fig.update_traces(marker=dict(opacity=plot_df['alpha']))

        pos_count = min(n_labels, len(sig_genes))
        if pos_count > 0:
            top_positive = sig_genes.nlargest(pos_count, 'R').index
            valid_indices = [idx for idx in top_positive if idx in plot_df.index]
            
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
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=40, b=40))

        return fig

    # --- TAB 2: SITE PROFILING (BAR CHART) ---
    @render_widget
    def site_profile_plot():
        target = input.target_site()
        if not target or target == "No Data Found" or df.empty:
            return px.bar(title="No Site Selected")

        # Get the row corresponding to the selected site
        row = df[df['Labels'] == target].iloc[0]
        
        # Extract the engagement R values across all available drugs
        drug_data = []
        for d in drugs:
            log_col = f'log2 {d} R'
            if log_col in df.columns:
                r_val = 2**row[log_col]
                drug_data.append({'Drug': d, 'R Value': r_val})
                
        bar_df = pd.DataFrame(drug_data).dropna()
        
        # Create a Bar Chart sorted by engagement
        bar_df = bar_df.sort_values('R Value', ascending=False)
        fig = px.bar(bar_df, x='Drug', y='R Value', title=f"Engagement Profile: {target}", color='R Value', color_continuous_scale='Reds')
        
        # Add the threshold line for reference
        fig.add_hline(y=input.threshold(), line_width=1, line_dash='dash', line_color='red', annotation_text="Threshold")
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=20, r=20, t=40, b=20), coloraxis_showscale=False)
        
        return fig

    # --- TAB 2: ALPHAFOLD VIEWER ---
    @render.ui
    def alphafold_viewer():
        target = input.target_site()
        if not target or target == "No Data Found" or df.empty:
            return ui.p("Please select a site from the dropdown.")

        # Extract UniProt ID and Site position
        row = df[df['Labels'] == target].iloc[0]
        uniprot_raw = str(row['Protein Id'])
        site_pos = str(row['Site Position'])
        
        # Clean the UniProt ID (AlphaFold URLs don't take isoforms like P04637-2)
        uniprot = uniprot_raw.split('-')[0]

        # Fetch PDB from AlphaFold Database
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb"
        try:
            req = urllib.request.urlopen(url)
            pdb_data = req.read().decode('utf-8')
        except Exception:
            return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve structure for {uniprot}. It may not exist in the AlphaFold database or the network request failed.</div>")

        # Initialize py3Dmol viewer
        view = py3Dmol.view(width="100%", height=500)
        view.addModel(pdb_data, "pdb")
        
        # Style the protein (rainbow scheme based on secondary structure/position)
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        
        # Highlight the specific site (stick + sphere)
        view.addStyle({'resi': site_pos}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
        view.addStyle({'resi': site_pos}, {'sphere': {'color': 'red', 'radius': 1.0}})
        
        # Zoom strictly to the highlighted site
        view.zoomTo({'resi': site_pos})
        
        # Embed securely
        viewer_html = html.escape(view._make_html())
        iframe = f'<iframe srcdoc="{viewer_html}" style="width: 100%; height: 500px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>'
        
        return ui.HTML(iframe)

# --- 4. Run App ---
app = App(app_ui, server)