import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
import html
import base64
import json
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

# --- NEW: Load the pre-calculated PPI database ---
ppi_db_path = app_dir / "ppi_reference.csv"
if ppi_db_path.exists():
    ppi_df = pd.read_csv(ppi_db_path)
else:
    ppi_df = pd.DataFrame(columns=['Target', 'PDB_ID', 'Min_Distance'])

# Safely load the Smiles dictionary & Cancer List
ypt_lib_path = app_dir / 'ypt_library.csv'
if ypt_lib_path.exists():
    ypt_lib = pd.read_csv(ypt_lib_path)
    smiles_dict = dict(zip(ypt_lib['Catalog ID'], ypt_lib['Smiles']))
    type_dict = dict(zip(ypt_lib['Catalog ID'], ypt_lib['Type']))
else:
    smiles_dict = {}
    type_dict = {}
    
cancer_path = app_dir / 'cancer_gene_shortlist.csv'
if cancer_path.exists():
    cancer = pd.read_csv(cancer_path)
else:
    cancer = pd.DataFrame(columns=['Gene']) 

if len(all_files) > 0:
    df_list = [pd.read_csv(f) for f in all_files]
    df = reduce(lambda left, right: pd.merge(left, right, on='Info', how='outer'), df_list)

    if 'Info' in df.columns:
        df[['Protein Id', 'Gene Symbol', 'Site Position', 'Sequence', 'Description']] = df['Info'].str.split('++', expand=True, regex=False)
        df.drop(columns='Info', inplace=True)
        df['Site Position'] = df['Site Position'].fillna('Unknown')
        df['Labels'] = df['Gene Symbol'] + "_Y" + df['Site Position'].astype(str)

    raw_drugs = list(set([i.split(' ')[1] for i in df.columns if 'log' in i]))
    
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
            
    sorted_drugs = sorted(raw_drugs, key=lambda x: drug_promiscuity[x])
    drug_choices = {d: f"{d} ({type_dict.get(d, 'Unknown')}, {drug_promiscuity[d]:.2f}%)" for d in sorted_drugs}
    default_drug = sorted_drugs[0] if sorted_drugs else None
    
    gene_to_sites = df.dropna(subset=['Gene Symbol', 'Site Position']).groupby('Gene Symbol')['Site Position'].apply(lambda x: sorted(list(set(x)))).to_dict()
    gene_choices = sorted(list(gene_to_sites.keys()))
    default_gene = gene_choices[0] if gene_choices else None
    default_site_choices = gene_to_sites.get(default_gene, [])
    default_site = default_site_choices[0] if default_site_choices else None
else:
    df = pd.DataFrame()
    drug_choices = {"No Data Found": "No Data Found"}
    default_drug = "No Data Found"
    gene_to_sites = {}
    gene_choices = ["No Data Found"]
    default_gene = "No Data Found"
    default_site_choices = ["No Data Found"]
    default_site = "No Data Found"

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Tyrosine Library Screening"),
    
    ui.navset_card_tab(
        ui.nav_panel(
            "Compound View", 
            ui.layout_columns(
                ui.div(
                    ui.card(
                        ui.h5("Settings"),
                        ui.input_selectize("data_type", "Select Drug (Type, Promiscuity):", choices=drug_choices, selected=default_drug),
                        ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
                        ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
                        ui.input_select(
                            "color_mode", "Color Points By:", 
                            choices=["Above Threshold", "P-Value Gradient", "Cancer-Driver List", "Highlight Custom List"]
                        ),
                        ui.input_text("custom_list", "Genes to Highlight (comma-separated):", placeholder="e.g. MAPK1, EGFR")
                    ),
                    ui.card(
                        ui.h5("Chemical Structure"),
                        ui.output_ui("molecule_ui_compound")
                    )
                ),
                ui.div(
                    ui.card(
                        ui.h5("Proteome Engagement Profile"),
                        output_widget("site_plot"),
                        ui.hr(),
                        ui.h5("Volcano Plot"),
                        output_widget("volcano_plot"),
                        style="padding: 10px; height: 100%; overflow-y: auto;"
                    )
                ),
                col_widths=(3, 9) 
            ),
            value="compound_tab"
        ),
        ui.nav_panel(
            "Site View",
            ui.layout_columns(
                ui.div(
                    ui.card(
                        ui.h5("Compound Engagement Profile"),
                        ui.layout_columns(
                            ui.input_selectize("target_gene", "Gene:", choices=gene_choices, selected=default_gene),
                            ui.input_selectize("target_site_pos", "Site:", choices=default_site_choices, selected=default_site),
                            col_widths=(6, 6)
                        ),
                        output_widget("site_profile_plot")
                    ),
                    ui.card(
                        ui.h5("Chemical Structure"),
                        ui.output_ui("molecule_ui_site")
                    )
                ),
                ui.div(
                    ui.navset_card_tab(
                        ui.nav_panel(
                            "Experimental Structure",
                            ui.div(
                                ui.p("Selected site is highlighted in red.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                                ui.input_select("pdb_selector", "Select PDB structure:", choices=["Loading..."], width="300px"),
                                style="display: flex; justify-content: space-between; align-items: center; padding-bottom: 10px;"
                            ),
                            ui.output_ui("pdb_viewer")
                        ),
                        ui.nav_panel(
                            "AlphaFold Structure",
                            ui.p("Selected site is highlighted in red.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                            ui.output_ui("alphafold_viewer")
                        ),
                        id="structure_tabs"
                    ),

                    # Inside app_ui, under the PPI card:
                    ui.card(
                        ui.h5("Protein-Protein Interactions (PPI)"),
                        ui.p("Ranked by proximity to site.", style="color: gray; font-size: 0.9em; margin-bottom: 10px;"),
                        ui.layout_columns(
                            ui.input_select("ppi_selector", "Select Interface (PDB | Distance):", choices=["Loading..."]),
                            col_widths=(12,) # Full width for the descriptive dropdown
                        ),
                        ui.output_ui("ppi_viewer")
                    )
                ),
                col_widths=(5, 7) 
            ),
            value="site_tab" 
        ),
        id="main_tabs"
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

    @reactive.Effect
    def update_site_dropdown():
        gene = input.target_gene()
        if gene and gene in gene_to_sites:
            sites = gene_to_sites[gene]
            ui.update_selectize("target_site_pos", choices=sites, selected=sites[0] if sites else None)

    @reactive.Calc
    def available_pdbs():
        gene = input.target_gene()
        site_str = input.target_site_pos()
        if not gene or not site_str or gene == "No Data Found" or df is None or df.empty:
            return {}

        target = f"{gene}_Y{site_str}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty:
            return {}

        row = matching_rows.iloc[0]
        protein_string = str(row['Protein Id'])
        uniprot = (protein_string.split('|')[1] if '|' in protein_string else protein_string).split('-')[0]

        try:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
            req = urllib.request.urlopen(uniprot_url)
            data = json.loads(req.read().decode('utf-8'))
            
            pdb_list = []
            site_num = int(site_str) if str(site_str).isdigit() else -1

            # 1. Parse and extract metadata
            for ref in data.get('uniProtKBCrossReferences', []):
                if ref['database'] == 'PDB':
                    pdb_id = ref['id']
                    method = "Unknown"
                    res_val = 999.0  # Default for sorting (worse than any real resolution)
                    has_site = False
                    coverage = 0
                    
                    for prop in ref.get('properties', []):
                        if prop['key'] == 'Method': 
                            method = prop['value']
                        elif prop['key'] == 'Resolution' and prop['value'] != '-':
                            try:
                                res_val = float(prop['value'].replace('A', '').strip())
                            except: pass
                        elif prop['key'] == 'Chains' and '=' in prop['value']:
                            ranges = prop['value'].split('=', 1)[1]
                            for r in ranges.split(','):
                                bounds = r.strip().split('-')
                                if len(bounds) >= 2 and bounds[0].isdigit() and bounds[-1].isdigit():
                                    start, end = int(bounds[0]), int(bounds[-1])
                                    coverage += (end - start + 1)
                                    if site_num != -1 and start <= site_num <= end:
                                        has_site = True
                    
                    # 2. Assign Method Score (Lower is better for sorting)
                    method_score = 4 # Default
                    if 'EM' in method.upper(): method_score = 1
                    elif 'X-RAY' in method.upper(): method_score = 2
                    elif 'NMR' in method.upper(): method_score = 3
                                        
                    pdb_list.append({
                        'id': pdb_id,
                        'has_site': has_site,
                        'method_score': method_score,
                        'res_val': res_val,
                        'coverage': coverage,
                        'method_label': method,
                        'res_label': "N/A" if res_val == 999.0 else f"{res_val}Å"
                    })
            
            # 3. Master Sort Logic
            # Sort keys: Site (True first), Method (Low score first), Resolution (Low first), Length (High first)
            pdb_list.sort(key=lambda x: (
                not x['has_site'],   # True (0) comes before False (1)
                x['method_score'],   # 1 (EM) comes before 2 (X-Ray)
                x['res_val'],        # 2.0 comes before 3.5
                -x['coverage']       # High coverage comes before low
            ))
            
            # 4. Format Labels for Dropdown
            pdb_choices = {}
            for p in pdb_list:
                site_tag = "Site Present" if p['has_site'] else "Site Absent"
                label = f"{p['id']}|{site_tag}|{p['method_label']}|{p['res_label']}|{p['coverage']}aa"
                pdb_choices[p['id']] = label
                
            return pdb_choices
        except Exception:
            return {}

    @reactive.Effect
    def update_pdb_dropdown():
        pdb_dict = available_pdbs()
        if pdb_dict:
            first_key = list(pdb_dict.keys())[0]
            ui.update_select("pdb_selector", choices=pdb_dict, selected=first_key)
        else:
            ui.update_select("pdb_selector", choices={"none": "No PDBs Found"}, selected="none")

    # --- RANKED PPI DROPDOWN LOGIC ---
    @reactive.Effect
    def update_ppi_dropdown():
        gene = input.target_gene()
        site_str = input.target_site_pos()

        if not gene or not site_str or ppi_df.empty:
            ui.update_select("ppi_selector", choices={"none": "No PPI Data Available"}, selected="none")
            return

        target = f"{gene}_Y{site_str}"
        
        # 1. Grab all structures involving this specific site
        valid_ppis = ppi_df[ppi_df['Target'] == target].copy()
        
        if valid_ppis.empty:
            ui.update_select("ppi_selector", choices={"none": "No structural interfaces found"}, selected="none")
        else:
            # 2. Sort strictly by distance (closest interactions at the top)
            valid_ppis = valid_ppis.sort_values('Min_Distance', ascending=True)
            
            # 3. Build descriptive labels including the distance
            choices = {}
            for _, row in valid_ppis.iterrows():
                dist = row['Min_Distance']
                label = f"{row['PDB_ID']}|{dist:.1f} Å"
                choices[row['PDB_ID']] = label
                
            first_key = list(choices.keys())[0]
            ui.update_select("ppi_selector", choices=choices, selected=first_key)

    # --- SHARED CHEMICAL HTML GENERATOR ---
    def generate_molecule_html(drug):
        if not drug:
            return ui.p("Select a drug or click a bar to view structure.", style="color: gray; font-style: italic;")
        smiles = smiles_dict.get(drug, "") 
        if not smiles:
            return ui.p("No SMILES available for this drug.", style="color: gray; font-style: italic;")
        safe_smiles = urllib.parse.quote(smiles)
        img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{safe_smiles}/PNG"
        return ui.HTML(f'<div style="text-align: center;"><img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Chemical Structure"></div>')

    @render.ui
    def molecule_ui_compound():
        return generate_molecule_html(active_drug.get())

    @render.ui
    def molecule_ui_site():
        return generate_molecule_html(active_drug.get())

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
            'R_plot': 2**df[f'log2 {data_type} R'],
            'R': 2**df[f'log2 {data_type} R'],
            'P-value': 10**-df[f'-log10 {data_type} p'],
            'log2 R': df[f'log2 {data_type} R'],
            '-log10 p': df[f'-log10 {data_type} p'],
            'Label': df['Labels'],
            'Gene Symbol': df['Gene Symbol'],
            'ID': df['Protein Id'],
            'Description': df['Description'],
        })
        
        plot_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        plot_df.dropna(subset=['R', 'P-value', 'log2 R', '-log10 p'], inplace=True)
        
        plot_df['Site Rank'] = plot_df['R'].rank(ascending=False).astype(int).astype(str)+'/'+str(len(plot_df))
        plot_df['% Site Rank'] = plot_df['R'].rank(ascending=False)/len(plot_df)*100

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

        valid_indices = []
        pos_count = min(n_labels, len(sig_genes))
        if pos_count > 0:
            top_positive = sig_genes.nlargest(pos_count, 'R').index
            valid_indices = [idx for idx in top_positive if idx in plot_df.index]
            
        return plot_df, valid_indices

    hover_dict = {
        'Site': False, 'color': False, 'R_plot': False,
        'ID': True, 'Label': True, 'Description': True,     
        'log2 R': False, '-log10 p': False, 'R': ':.3f', 'P-value': ':.4f',
        'Site Rank': True, '% Site Rank': ':.5f'
    }

    # --- PROFILE PLOT ---
    @render_widget
    def site_plot():
        plot_df, valid_indices = plot_data()
        color_mode = input.color_mode()
        threshold = input.threshold()

        if color_mode == "Above Threshold":
            fig = px.scatter(
                plot_df, x='Site', y='R_plot', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=0.9), selector=dict(name='high'))
        elif color_mode in ["Highlight Custom List", "Cancer-Driver List"]:
            fig = px.scatter(
                plot_df, x='Site', y='R_plot', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'highlight': '#D62728'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.15), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))
        else: # P-Value Gradient
            fig = px.scatter(
                plot_df, x='Site', y='R_plot', color='P-value',
                color_continuous_scale='Viridis_r', 
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=plot_df['alpha']))

        for idx in valid_indices:
            row = plot_df.loc[idx]
            fig.add_annotation(
                x=row['Site'], y=row['R'], text=row['Label'],
                showarrow=True, arrowsize=1, arrowwidth=1,
                arrowcolor="#444", ax=20, ay=-30,
                font=dict(size=10, color="black"), bgcolor="white",
                bordercolor="#c7c7c7", borderwidth=1, borderpad=3
            )

        fig.add_hline(y=threshold, line_width=1, line_dash='dash')
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(yaxis_title="R", plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))

        widget = go.FigureWidget(fig)
        current_config = widget._config or {}
        widget._config = {**current_config, 'edits': {'annotationTail': True}}
        return widget

    # --- VOLCANO PLOT ---
    @render_widget
    def volcano_plot():
        plot_df, _ = plot_data() 
        volcano_df = plot_df.copy()
        color_mode = input.color_mode()
        threshold = input.threshold()

        safe_thresh = threshold if threshold > 0 else 1e-9 
        sig_mask = (volcano_df['log2 R'] > np.log2(safe_thresh)) & (volcano_df['-log10 p'] > -np.log10(0.05))

        volcano_df['color'] = 'non-significant'

        if color_mode == "Highlight Custom List":
            custom_genes = [g.strip().upper() for g in input.custom_list().split(',') if g.strip()]
            mask = volcano_df['Gene Symbol'].str.upper().isin(custom_genes)
            volcano_df.loc[sig_mask & mask, 'color'] = 'highlight'
        elif color_mode == "Cancer-Driver List":
            mask = volcano_df['Gene Symbol'].isin(cancer['Gene'])
            volcano_df.loc[sig_mask & mask, 'color'] = 'highlight'
        else: 
            volcano_df.loc[sig_mask, 'color'] = 'high'

        colored_genes = volcano_df[volcano_df['color'] != 'non-significant']
        pos_count = min(input.n_labels(), len(colored_genes))
        
        valid_indices = []
        if pos_count > 0:
            top_positive = colored_genes.nlargest(pos_count, 'R').index
            valid_indices = top_positive.tolist()

        if color_mode == "P-Value Gradient":
            fig = px.scatter(
                volcano_df[volcano_df['color'] == 'non-significant'], 
                x='log2 R', y='-log10 p', color_discrete_sequence=['#dddddd'], hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.2), hovertemplate=None)
            
            sig_df = volcano_df[volcano_df['color'] == 'high']
            if not sig_df.empty:
                fig2 = px.scatter(
                    sig_df, x='log2 R', y='-log10 p', color='P-value', 
                    color_continuous_scale='Viridis_r', hover_data=hover_dict
                )
                for trace in fig2.data: fig.add_trace(trace)
                fig.layout.coloraxis = fig2.layout.coloraxis
        else:
            fig = px.scatter(
                volcano_df, x='log2 R', y='-log10 p', color='color',
                color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD', 'highlight': '#D62728'},
                hover_data=hover_dict
            )
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))
            fig.update_traces(marker=dict(opacity=0.9), selector=dict(name='high'))
            fig.update_traces(marker=dict(opacity=1.0), selector=dict(name='highlight'))

        for idx in valid_indices:
            row = volcano_df.loc[idx]
            fig.add_annotation(
                x=row['log2 R'], y=row['-log10 p'], text=row['Label'],
                showarrow=True, arrowsize=1, arrowwidth=1, arrowcolor="#444", ax=20, ay=-30,
                font=dict(size=10, color="black"), bgcolor="white", bordercolor="#c7c7c7", borderwidth=1, borderpad=3
            )

        fig.add_vline(x=np.log2(safe_thresh), line_width=1, line_dash='dash')
        fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash='dash')
        fig.update_traces(marker=dict(size=6))
        fig.update_layout(xaxis_title="log<sub>2</sub> Fold Change", yaxis_title="-log<sub>10</sub> P-Value", plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))

        widget = go.FigureWidget(fig)
        current_config = widget._config or {}
        widget._config = {**current_config, 'edits': {'annotationTail': True}}
        return widget

    # --- SITE PROFILING ---
    @render_widget
    def site_profile_plot():
        gene = input.target_gene()
        site = input.target_site_pos()
        
        if not gene or not site or gene == "No Data Found" or df.empty:
            return go.FigureWidget(px.bar(title="No Site Selected"))

        target = f"{gene}_Y{site}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty: return go.FigureWidget(px.bar(title="Loading..."))
             
        row = matching_rows.iloc[0]
        
        drug_data = []
        for d in raw_drugs: 
            log_col = f'log2 {d} R'
            if log_col in df.columns:
                r_val = 2**row[log_col]
                drug_data.append({'Drug': d, 'R': r_val})
                
        bar_df = pd.DataFrame(drug_data).dropna().sort_values('R', ascending=False)
        fig = px.bar(bar_df, x='Drug', y='R', color='R', color_continuous_scale='Reds', hover_data = {'Drug': True, 'R': ':.3f'})
        fig.add_hline(y=input.threshold(), line_width=1, line_dash='dash', line_color='red', annotation_text="Threshold")
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=20, r=20, t=40, b=20), coloraxis_showscale=False)
        
        widget = go.FigureWidget(fig)
        
        def handle_bar_click(trace, points, state):
            if points.point_inds:
                idx = points.point_inds[0]
                clicked_drug = trace.x[idx]
                active_drug.set(clicked_drug)
                
        if widget.data: widget.data[0].on_click(handle_bar_click)
        return widget

    # --- ALPHAFOLD VIEWER ---
    @render.ui
    def alphafold_viewer():
        gene = input.target_gene()
        site = input.target_site_pos()
        if not gene or not site or gene == "No Data Found" or df.empty: return ui.p("Please select a Gene and Site.")

        target = f"{gene}_Y{site}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty: return ui.p("Loading structure...")

        row = matching_rows.iloc[0]
        protein_string = str(row['Protein Id'])
        uniprot = (protein_string.split('|')[1] if '|' in protein_string else protein_string).split('-')[0]
        site_pos = str(row['Site Position'])

        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.cif"
        try:
            pdb_data = urllib.request.urlopen(url).read().decode('utf-8')
            pdb_data = pdb_data.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception:
            return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve structure for {uniprot}.</div>")

        view = py3Dmol.view(width="100%", height=500)
        view.addModel(pdb_data, "cif")
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.addStyle(
            {'resi': site_pos, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, 
            {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}}
        )
        view.zoomTo({'resi': site_pos})

        b64_html = base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')
        return ui.HTML(f'<iframe src="data:text/html;base64,{b64_html}" style="width: 100%; height: 500px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>')
    
    # --- PDB VIEWER ---
    @render.ui
    def pdb_viewer():
        pdb_id = input.pdb_selector()
        if not pdb_id or pdb_id in ["Loading...", "none"]: return ui.HTML("<div style='color:orange; padding:20px;'><b>Note:</b> No PDB structure found.</div>")
        if pdb_id not in available_pdbs(): return ui.p("Syncing...")

        gene = input.target_gene()
        site = input.target_site_pos()
        target = f"{gene}_Y{site}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty: return ui.p("")
        site_pos = str(matching_rows.iloc[0]['Site Position'])

        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        try:
            pdb_data = urllib.request.urlopen(url).read().decode('utf-8')
            pdb_data = pdb_data.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception:
            return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve PDB {pdb_id}.</div>")

        view = py3Dmol.view(width="100%", height=480)
        view.addModel(pdb_data, "cif")
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.addStyle(
            {'resi': site_pos, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, 
            {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}}
        )
        view.zoomTo({'resi': site_pos})
        
        b64_html = base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')
        return ui.HTML(f'''
        <div style="margin-bottom: 5px; font-size: 0.9em; line-height: 1.3;">
            <b>Displaying PDB: <a href="https://www.rcsb.org/structure/{pdb_id}" target="_blank">{pdb_id}</a></b>
        </div>
        <iframe src="data:text/html;base64,{b64_html}" style="width: 100%; height: 440px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>
        ''')

    # --- PPI VIEWER ---
    @render.ui
    def ppi_viewer():
        pdb_id = input.ppi_selector()
        target_site = input.target_site_pos()
        
        if not pdb_id or pdb_id in ["none", "Loading..."]:
            return ui.div()

        # Retrieve the distance from the DB to set the viewer scale
        target_label = f"{input.target_gene()}_Y{target_site}"
        match = ppi_df[(ppi_df['Target'] == target_label) & (ppi_df['PDB_ID'] == pdb_id)]
        
        # Default to 5.0 if not found, otherwise add a 1A buffer to the visual highlight
        viz_cutoff = float(match['Min_Distance'].iloc[0]) + 1.0 if not match.empty else 5.0

        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        try:
            pdb_data = urllib.request.urlopen(url).read().decode('utf-8')
            # Essential sanitization for JS safety
            pdb_data = pdb_data.replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception as e:
            return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve {pdb_id}.</div>")

        view = py3Dmol.view(width="100%", height=480)
        view.addModel(pdb_data, "cif")
        
        # Biological coloring
        view.setStyle({'cartoon': {'colorscheme': 'chain'}})

        # Sidechain-only highlight for the target
        view.addStyle(
            {'resi': target_site, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, 
            {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}}
        )

        # Highlight the partner chain atoms that are within range
        view.addStyle(
            {'within': {'distance': viz_cutoff, 'sel': {'resi': target_site}}, 'not': {'atom': ['N', 'C', 'O', 'OXT']}},
            {'stick': {'colorscheme': 'cyanCarbon'}}
        )
        
        view.zoomTo({'resi': target_site})
        
        b64_html = base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')

        return ui.HTML(f'''
        <div style="margin-bottom: 5px; font-size: 0.9em; line-height: 1.3;">
            <b>PDB: <a href="https://www.rcsb.org/structure/{pdb_id}" target="_blank">{pdb_id}</a></b><br>
            <span style="color: #444;"><i>Interface detected at {viz_cutoff-1.0:.1f} Å. Partner chains shown in Cyan.</i></span>
        </div>
        <iframe src="data:text/html;base64,{b64_html}" style="width: 100%; height: 440px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>
        ''')

# --- 4. Run App ---
app = App(app_ui, server)