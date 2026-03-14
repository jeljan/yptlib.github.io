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

# Load PPI database
ppi_db_path = app_dir / "ppi_reference.csv"
if ppi_db_path.exists():
    ppi_df = pd.read_csv(ppi_db_path)
else:
    ppi_df = pd.DataFrame(columns=['Target', 'PDB_ID', 'Min_Distance'])

# Load Smiles dictionary & Cancer List
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
    drug_choices = {d: f"{d} ({type_dict.get(d, 'Unknown')}, {drug_promiscuity[d]:.2f}%)"} # truncated for display
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
                        ui.input_selectize("data_type", "Select Drug:", choices=drug_choices, selected=default_drug),
                        ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
                        ui.input_numeric("n_labels", "Number of Labels:", value=5, min=0, max=20),
                        ui.input_select("color_mode", "Color Points By:", choices=["Above Threshold", "P-Value Gradient", "Cancer-Driver List", "Highlight Custom List"]),
                        ui.input_text("custom_list", "Genes to Highlight:", placeholder="e.g. MAPK1, EGFR")
                    ),
                    ui.card(ui.h5("Chemical Structure"), ui.output_ui("molecule_ui_compound"))
                ),
                ui.div(
                    ui.card(ui.h5("Proteome Engagement Profile"), output_widget("site_plot"), ui.hr(), ui.h5("Volcano Plot"), output_widget("volcano_plot"))
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
                    ui.card(ui.h5("Chemical Structure"), ui.output_ui("molecule_ui_site"))
                ),
                ui.div(
                    # --- MASTER STRUCTURAL CONTROL PANEL ---
                    ui.div(
                        ui.input_switch("spin_toggle", "Auto-Spin", value=False),
                        ui.input_select("style_choice", None, choices={"cartoon": "Cartoon Style", "sphere": "Sphere Style", "line": "Line Style"}, width="160px"),
                        ui.input_action_button("zoom_btn", "Recenter View", class_="btn-sm btn-outline-secondary"),
                        style="display: flex; gap: 20px; align-items: center; padding: 10px; background: #fdfdfd; border: 1px solid #eee; border-radius: 8px 8px 0 0; margin-bottom: -1px;"
                    ),
                    ui.navset_card_tab(
                        ui.nav_panel(
                            "Experimental Structure",
                            ui.div(
                                ui.input_select("pdb_selector", None, choices=["Loading..."], width="350px"),
                                style="padding: 10px; border-bottom: 1px solid #eee;"
                            ),
                            ui.output_ui("pdb_viewer")
                        ),
                        ui.nav_panel(
                            "AlphaFold Structure",
                            ui.output_ui("alphafold_viewer")
                        ),
                        ui.nav_panel(
                            "PPI Interfaces",
                            ui.div(
                                ui.input_select("ppi_selector", None, choices=["Loading..."], width="350px"),
                                style="padding: 10px; border-bottom: 1px solid #eee;"
                            ),
                            ui.output_ui("ppi_viewer")
                        ),
                        id="structure_tabs"
                    ),
                    # --- SEQUENCE NAVIGATOR ---
                    ui.card(
                        ui.h6("Protein Sequence Navigator", style="font-size: 0.8rem; color: #666; margin-bottom: 5px;"),
                        ui.div(
                            ui.output_ui("sequence_bar"),
                            style="overflow-x: auto; white-space: nowrap; padding: 12px; background: #fff; border: 1px solid #eee; border-radius: 4px; font-family: 'Courier New', Courier, monospace; font-size: 14px;"
                        )
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
        if input.main_tabs() == "site_tab": active_drug.set(None) 
        else: active_drug.set(input.data_type()) 

    @reactive.Effect
    def update_site_dropdown():
        gene = input.target_gene()
        if gene and gene in gene_to_sites:
            sites = gene_to_sites[gene]
            ui.update_selectize("target_site_pos", choices=sites, selected=sites[0] if sites else None)

    # --- PDB RANKING LOGIC ---
    @reactive.Calc
    def available_pdbs():
        gene, site_str = input.target_gene(), input.target_site_pos()
        if not gene or not site_str or gene == "No Data Found" or df is None: return {}
        target = f"{gene}_Y{site_str}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty: return {}
        row = matching_rows.iloc[0]
        uniprot = str(row['Protein Id']).split('|')[1].split('-')[0] if '|' in str(row['Protein Id']) else str(row['Protein Id']).split('-')[0]

        try:
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
            data = json.loads(urllib.request.urlopen(uniprot_url).read().decode('utf-8'))
            pdb_list = []
            site_num = int(site_str) if site_str.isdigit() else -1

            for ref in data.get('uniProtKBCrossReferences', []):
                if ref['database'] == 'PDB':
                    pdb_id = ref['id']
                    method, res_val, has_site, coverage = "Unknown", 999.0, False, 0
                    for prop in ref.get('properties', []):
                        if prop['key'] == 'Method': method = prop['value']
                        elif prop['key'] == 'Resolution':
                            try: res_val = float(prop['value'].replace('A', ''))
                            except: pass
                        elif prop['key'] == 'Chains' and '=' in prop['value']:
                            ranges = prop['value'].split('=', 1)[1]
                            for r in ranges.split(','):
                                bounds = r.strip().split('-')
                                if len(bounds) >= 2:
                                    start, end = int(bounds[0]), int(bounds[-1])
                                    coverage += (end - start + 1)
                                    if start <= site_num <= end: has_site = True
                    
                    m_score = 1 if 'EM' in method.upper() else 2 if 'X-RAY' in method.upper() else 3 if 'NMR' in method.upper() else 4
                    pdb_list.append({'id': pdb_id, 'has_site': has_site, 'm_score': m_score, 'res': res_val, 'cov': coverage, 'm_label': method})
            
            pdb_list.sort(key=lambda x: (not x['has_site'], x['m_score'], x['res'], -x['cov']))
            return {p['id']: f"{'★' if p['has_site'] else '○'} {p['id']} | {p['m_label']} | {p['res']}Å | {p['cov']}aa" for p in pdb_list}
        except: return {}

    @reactive.Effect
    def update_pdb_dropdown():
        pdb_dict = available_pdbs()
        if pdb_dict: ui.update_select("pdb_selector", choices=pdb_dict, selected=list(pdb_dict.keys())[0])
        else: ui.update_select("pdb_selector", choices={"none": "No PDBs Found"}, selected="none")

    @reactive.Effect
    def update_ppi_dropdown():
        target = f"{input.target_gene()}_Y{input.target_site_pos()}"
        valid_ppis = ppi_df[ppi_df['Target'] == target].sort_values('Min_Distance')
        if valid_ppis.empty: ui.update_select("ppi_selector", choices={"none": "No interfaces found"}, selected="none")
        else:
            choices = {row['PDB_ID']: f"{'⚡' if row['Min_Distance'] < 5 else '🟢'} {row['PDB_ID']} | {row['Min_Distance']:.1f} Å" for _, row in valid_ppis.iterrows()}
            ui.update_select("ppi_selector", choices=choices, selected=list(choices.keys())[0])

    # --- SEQUENCE BAR LOGIC ---
    @render.ui
    def sequence_bar():
        target = f"{input.target_gene()}_Y{input.target_site_pos()}"
        matching = df[df['Labels'] == target]
        if matching.empty: return ui.p("No sequence.")
        
        full_seq = str(matching.iloc[0]['Sequence'])
        site = int(input.target_site_pos())
        seq_html = []
        for i, aa in enumerate(full_seq):
            res_num = i + 1
            is_hit = res_num == site
            style = f"color:{'#d62728' if is_hit else '#555'}; font-weight:{'bold' if is_hit else 'normal'}; border:{'1px solid #d62728' if is_hit else 'none'}; padding:2px; display:inline-block; min-width:12px; text-align:center;"
            seq_html.append(f'<span title="Residue {res_num}" style="{style}">{aa}</span>')
        return ui.HTML("".join(seq_html))

    # --- SHARED 3D RENDER ENGINE ---
    def render_structure(pdb_id, is_alphafold=False, is_ppi=False):
        site_pos = input.target_site_pos()
        if is_alphafold:
            uniprot = str(df[df['Labels'] == f"{input.target_gene()}_Y{site_pos}"].iloc[0]['Protein Id']).split('|')[1].split('-')[0]
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.cif"
        else:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
        
        try:
            pdb_data = urllib.request.urlopen(url).read().decode('utf-8').replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
            view = py3Dmol.view(width="100%", height=480)
            view.addModel(pdb_data, "cif")
            
            # Apply Style Control
            style = input.style_choice()
            view.setStyle({style: {'colorscheme': 'chain' if is_ppi else 'spectrum'}})
            
            # Sidechain-only Highlight (Excluding Backbone)
            view.addStyle(
                {'resi': site_pos, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, 
                {'stick': {'colorscheme': 'redCarbon', 'radius': 0.3}}
            )

            if is_ppi:
                # Highlight interface atoms in Cyan
                match = ppi_df[(ppi_df['Target'] == f"{input.target_gene()}_Y{site_pos}") & (ppi_df['PDB_ID'] == pdb_id)]
                dist = float(match['Min_Distance'].iloc[0]) + 1.0 if not match.empty else 6.0
                view.addStyle(
                    {'within': {'distance': dist, 'sel': {'resi': site_pos}}, 'not': {'atom': ['N', 'C', 'O', 'OXT']}},
                    {'stick': {'colorscheme': 'cyanCarbon'}}
                )

            view.zoomTo({'resi': site_pos})
            if input.spin_toggle(): view.spin(True)
            
            b64 = base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')
            return ui.HTML(f'<iframe src="data:text/html;base64,{b64}" style="width:100%; height:480px; border:none; overflow:hidden;"></iframe>')
        except: return ui.p("Error loading structure.")

    @render.ui
    def alphafold_viewer(): return render_structure(None, is_alphafold=True)

    @render.ui
    def pdb_viewer():
        pid = input.pdb_selector()
        return render_structure(pid) if pid and pid != "none" else ui.div()

    @render.ui
    def ppi_viewer():
        pid = input.ppi_selector()
        return render_structure(pid, is_ppi=True) if pid and pid != "none" else ui.div()

    # --- Plotting logic (preserved from your snippet) ---
    @reactive.Calc
    def plot_data():
        dt, th, nl, cm = input.data_type(), input.threshold(), input.n_labels(), input.color_mode()
        p_df = pd.DataFrame({'Site': np.arange(len(df[f'log2 {dt} R'])), 'R': 2**df[f'log2 {dt} R'], 'P': 10**-df[f'-log10 {dt} p'], 'Label': df['Labels'], 'Gene': df['Gene Symbol']}).dropna()
        p_df['color'] = 'non-significant'
        if cm == "Above Threshold": p_df.loc[p_df['R'] > th, 'color'] = 'high'
        return p_df, p_df[p_df['R'] > th].nlargest(nl, 'R').index

    @render_widget
    def site_plot():
        d, ids = plot_data()
        fig = px.scatter(d, x='Site', y='R', color='color', color_discrete_map={'non-significant': '#ddd', 'high': '#4470AD'})
        return go.FigureWidget(fig)

    @render_widget
    def volcano_plot():
        d, _ = plot_data()
        fig = px.scatter(d, x=np.log2(d['R']), y=-np.log10(d['P']), color='color')
        return go.FigureWidget(fig)

    @render_widget
    def site_profile_plot():
        matching = df[df['Labels'] == f"{input.target_gene()}_Y{input.target_site_pos()}"]
        if matching.empty: return go.FigureWidget()
        row = matching.iloc[0]
        bar_data = [{'Drug': dr, 'R': 2**row[f'log2 {dr} R']} for dr in raw_drugs if f'log2 {dr} R' in df.columns]
        fig = px.bar(pd.DataFrame(bar_data).sort_values('R', ascending=False), x='Drug', y='R')
        return go.FigureWidget(fig)

    def generate_molecule_html(drug):
        if not drug or drug not in smiles_dict: return ui.p("No structure.")
        img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles_dict[drug])}/PNG"
        return ui.HTML(f'<img src="{img_url}" style="width:100%; background:white; border-radius:4px;">')

    @render.ui
    def molecule_ui_compound(): return generate_molecule_html(active_drug.get())
    
    @render.ui
    def molecule_ui_site(): return generate_molecule_html(active_drug.get())

app = App(app_ui, server)