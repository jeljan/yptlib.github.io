import numpy as np
import pandas as pd
import urllib.parse
import urllib.request
import html
import base64
import json
import py3Dmol
import re
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

# --- Load the pre-calculated PPI database ---
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
        df['Sequence'] = df['Sequence'].str[2:-2].replace(r'\W', '', regex=True)
        df.drop(columns='Info', inplace=True)
        df['Site Position'] = df['Site Position'].fillna('Unknown')
        df['Labels'] = df['Gene Symbol'] + "_Y" + df['Site Position'].astype(str)

    raw_drugs = list(set([i.split(' ')[1] for i in df.columns if 'log' in i]))
    
    drug_promiscuity = {}
    for d in raw_drugs:
        log_col = f'log2 {d} R'
        if log_col in df.columns:
            total_valid_sites = df[log_col].notna().sum()
            promiscuity = (df[log_col] > 1).sum() / total_valid_sites * 100 if total_valid_sites > 0 else 0.0
            drug_promiscuity[d] = promiscuity
            
    sorted_drugs = sorted(raw_drugs, key=lambda x: drug_promiscuity[x])
    drug_choices = {d: f"{d} ({type_dict.get(d, 'Unknown')}, {drug_promiscuity[d]:.2f}%)" for d in sorted_drugs}
    default_drug = sorted_drugs[0] if sorted_drugs else None
    
    def site_sort_key(site_str):
        match = re.search(r'\d+', str(site_str))
        return (int(match.group()), str(site_str)) if match else (float('inf'), str(site_str))
    
    gene_to_sites = df.dropna(subset=['Gene Symbol', 'Site Position']).groupby('Gene Symbol')['Site Position'].apply(lambda x: sorted(list(set(x)), key=site_sort_key)).to_dict()
    gene_choices = sorted(list(gene_to_sites.keys()))
    default_gene = gene_choices[0] if gene_choices else None
    default_site_choices = gene_to_sites.get(default_gene, [])
    default_site = default_site_choices[0] if default_site_choices else None
else:
    df = pd.DataFrame()
    raw_drugs = []
    drug_choices = {"No Data Found": "No Data Found"}
    default_drug = "No Data Found"
    gene_to_sites, gene_choices, default_gene, default_site_choices, default_site = {}, ["No Data Found"], "No Data Found", ["No Data Found"], "No Data Found"

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Tyrosine Library Screening"),
    ui.navset_card_tab(
        ui.nav_panel(
            "Summary View",
            ui.layout_columns(
                ui.card(ui.h5("Site Promiscuity"), output_widget("summary_site_hist")),
                ui.card(ui.h5("Compound Promiscuity by Type"), output_widget("summary_drug_hist")),
                col_widths=(6, 6)
            ),
            ui.card(
                ui.h5("Target Engagement"),
                ui.div(
                    ui.input_switch("summary_sig_only", "Significant only (p < 0.05)", value=False)
                ),
                ui.layout_columns(
                    ui.div(
                        ui.input_selectize("summary_cancer_site", "Cancer Driver Sites:", choices=["Loading..."]),
                        output_widget("summary_cancer_bar")
                    ),
                    ui.div(
                        ui.input_selectize("summary_ppi_site", "PPI Interface Sites:", choices=["Loading..."]),
                        output_widget("summary_ppi_bar")
                    ),
                    col_widths=(6, 6)
                ),
                ui.div(
                    ui.h6("Compound Structure", style="text-align: center; color: #555;"),
                    ui.output_ui("molecule_ui_summary"),
                    style="display: flex; flex-direction: column; align-items: center; justify-content: center; padding-top: 10px;"
                )
            ),
            value="summary_tab"
        ),
        ui.nav_panel(
            "Compound View", 
            ui.layout_columns(
                ui.div(
                    ui.card(
                        ui.h5("Settings"),
                        ui.input_selectize("data_type", "Select Drug (Type, Promiscuity):", choices=drug_choices, selected=default_drug),
                        ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
                        ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
                        ui.input_select("color_mode", "Color Points By:", choices=["Above Threshold", "P-Value Gradient", "Cancer-Driver List", "Sites with PPIs", "Highlight Custom List"]),
                        ui.input_text("custom_list", "Genes to Highlight (comma-separated):", placeholder="e.g. MAPK1, EGFR")
                    ),
                    ui.card(ui.h5("Compound Structure"), ui.output_ui("molecule_ui_compound"))
                ),
                ui.div(
                    ui.card(
                        ui.h5("Proteome Engagement Profile"), output_widget("site_plot"), ui.hr(),
                        ui.h5("Volcano Plot"), output_widget("volcano_plot"), style="padding: 10px; height: 100%; overflow-y: auto;"
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
                    ui.card(ui.h5("Compound Structure"), ui.output_ui("molecule_ui_site"))
                ),
                ui.div(
                    ui.navset_card_tab(
                        ui.nav_panel(
                            "Experimental Structure",
                            ui.div(
                                ui.p("Selected site highlighted in red, hover over residue to see more info.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                                ui.input_select("pdb_selector", "Select PDB structure:", choices=["Loading..."], width="350px"),
                                style="display: flex; justify-content: space-between; align-items: center; padding-bottom: 10px;"
                            ),
                            ui.output_ui("pdb_viewer")
                        ),
                        ui.nav_panel(
                            "AlphaFold Structure",
                            ui.p("Selected site highlighted in red, hover over residue to see more info.", style="color: gray; font-size: 0.9em; margin-bottom: 0;"),
                            ui.output_ui("alphafold_viewer")
                        ),
                        id="structure_tabs"
                    ),
                    ui.card(
                        ui.h5("Protein-Protein Interaction (PPI) Interfaces"),
                        ui.p("Selected site highlighted in red, hover over residue to see more info.", style="color: gray; font-size: 0.9em; margin-bottom: 10px;"),
                        ui.input_select("ppi_selector", "Select Interface (ID|Distance):", choices=["Loading..."]),
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

# --- MAPPING HELPERS ---
def get_chain_names(pdb_id):
    try:
        req = urllib.request.Request(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_id.lower()}")
        data = json.loads(urllib.request.urlopen(req).read().decode('utf-8'))
        chain_map = {}
        for mol in data.get(pdb_id.lower(), []):
            name = mol.get('molecule_name', ['Unknown Molecule'])[0]
            mol_type = mol.get('molecule_type', '').lower()
            if 'water' in mol_type or 'bound' in mol_type: continue
            for chain in mol.get('in_chains', []):
                if chain not in chain_map: chain_map[chain] = name
        return chain_map
    except Exception: return {}

def verify_and_map_site(pdb_id, site_pos, user_seq, pdbe_data=None):
    y_idx = str(user_seq).find('Y')
    if y_idx == -1: return None, None 

    if pdbe_data is None:
        try:
            req = urllib.request.Request(f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/{pdb_id.lower()}")
            pdbe_data = json.loads(urllib.request.urlopen(req, timeout=10).read().decode('utf-8'))
        except Exception:
            return 'A', str(site_pos)
            
    aa_map = {'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D', 'CYS':'C', 'GLN':'Q', 'GLU':'E', 
              'GLY':'G', 'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K', 'MET':'M', 'PHE':'F', 
              'PRO':'P', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V', 
              'PTR':'Y', 'TYS':'Y'}
              
    all_chains = []
    for mol in pdbe_data.get(pdb_id.lower(), {}).get('molecules', []):
        for chain in mol.get('chains', []):
            chain_id = chain.get('chain_id')
            chain_seq, auth_nums, observed = "", [], []
            for res in chain.get('residues', []):
                chain_seq += aa_map.get(res.get('residue_name', ''), 'X')
                auth_num = res.get('author_residue_number')
                auth_ins = res.get('author_insertion_code', '')
                auth_nums.append(f"{auth_num}{auth_ins}".strip() if auth_num is not None else None)
                observed.append(res.get('observed_ratio', 1) > 0)
                
            if chain_seq: all_chains.append((chain_id, chain_seq, auth_nums, observed))
            
    for chain_id, chain_seq, auth_nums, observed_flags in all_chains:
        start = 0
        while True:
            match_start = chain_seq.find(user_seq, start)
            if match_start == -1: break
            
            exact_match_idx = match_start + y_idx
            if exact_match_idx < len(auth_nums) and auth_nums[exact_match_idx] is not None:
                if observed_flags[exact_match_idx]: 
                    return chain_id, auth_nums[exact_match_idx]
            start = match_start + 1

    for chain_id, chain_seq, auth_nums, observed_flags in all_chains:
        if str(site_pos) in auth_nums:
            idx = auth_nums.index(str(site_pos))
            if idx < len(chain_seq) and chain_seq[idx] == 'Y':
                seq_start = max(0, idx - y_idx)
                seq_end = min(len(chain_seq), idx + len(user_seq) - y_idx)
                extracted = chain_seq[seq_start:seq_end]
                
                user_start = y_idx - (idx - seq_start)
                user_end = y_idx + (seq_end - idx)
                expected_subseq = user_seq[user_start:user_end]
                
                if extracted == expected_subseq:
                    if observed_flags[idx]: 
                        return chain_id, str(site_pos)

    return None, None

# --- 3. Shiny Server ---
def server(input, output, session):
    
    active_drug = reactive.Value(default_drug)
    
    unhover_js = """function(atom,viewer) {
        if(atom.label) {
            viewer.removeLabel(atom.label);
            delete atom.label;
        }
    }"""
        
    @reactive.Effect
    @reactive.event(input.main_tabs)
    def handle_tab_switch():
        if input.main_tabs() == "site_tab": active_drug.set(None)
        elif input.main_tabs() == "compound_tab": active_drug.set(input.data_type())
        elif input.main_tabs() == "summary_tab": active_drug.set(None)

    @reactive.Effect
    @reactive.event(input.data_type, ignore_init=True)
    def update_drug_from_dropdown():
        if input.main_tabs() == "compound_tab":
            active_drug.set(input.data_type())

    @reactive.Effect
    def update_site_dropdown():
        gene = input.target_gene()
        if gene and gene in gene_to_sites:
            sites = gene_to_sites[gene]
            ui.update_selectize("target_site_pos", choices=sites, selected=sites[0] if sites else None)

    # --- SUMMARY VIEW LOGIC ---
    @reactive.Calc
    def site_max_r():
        sig_only = input.summary_sig_only()
        if df is None or df.empty or not raw_drugs:
            return {}, {}
            
        r_cols = [f'log2 {d} R' for d in raw_drugs if f'log2 {d} R' in df.columns]
        r_df = df[r_cols].copy()
        
        if sig_only:
            p_thresh = -np.log10(0.05)
            for d in raw_drugs:
                r_col = f'log2 {d} R'
                p_col = f'-log10 {d} p'
                if r_col in df.columns and p_col in df.columns:
                    mask = (df[p_col].isna()) | (df[p_col] <= p_thresh)
                    r_df.loc[mask, r_col] = np.nan
                    
        r_df = 2 ** r_df
        max_r = r_df.max(axis=1)
        
        sum_df = pd.DataFrame({
            'Labels': df['Labels'],
            'Gene': df['Gene Symbol'],
            'Max_R': max_r
        }).dropna(subset=['Max_R'])
        
        cancer_df = sum_df[sum_df['Gene'].isin(cancer['Gene'])].sort_values('Max_R', ascending=False)
        cancer_choices = {row['Labels']: f"{row['Labels']} (Max R: {row['Max_R']:.2f})" for _, row in cancer_df.iterrows()}

        ppi_df_filtered = sum_df[sum_df['Labels'].isin(ppi_df['Target'])].sort_values('Max_R', ascending=False)
        ppi_choices = {row['Labels']: f"{row['Labels']} (Max R: {row['Max_R']:.2f})" for _, row in ppi_df_filtered.iterrows()}
        
        return cancer_choices, ppi_choices
  
    @reactive.Effect
    def update_summary_dropdowns():
        cancer_choices, ppi_choices = site_max_r()

        ui.update_selectize("summary_cancer_site", choices=cancer_choices, selected=list(cancer_choices.keys())[0] if cancer_choices else None)
        ui.update_selectize("summary_ppi_site", choices=ppi_choices, selected=list(ppi_choices.keys())[0] if ppi_choices else None)

    @render_widget
    def summary_site_hist():
        if df is None or df.empty: return go.FigureWidget()
        log_cols = [c for c in df.columns if 'log2' in c and ' R' in c]
        if not log_cols: return go.FigureWidget()
        
        site_prom = (df[[i for i in df.columns if 'R' in i]]>1).sum(axis=1)/df[[i for i in df.columns if 'R' in i]].shape[1]
        
        prom_df = pd.DataFrame({'Promiscuity': site_prom})
        fig = px.histogram(prom_df, x='Promiscuity', nbins=100, opacity=0.7, 
                           labels={'Promiscuity': 'Site Promiscuity (% Compounds Hit at R > 2)', 'Frequency': 'Frequency'})
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=40, r=40, t=10, b=40), legend_title_text='')
        return go.FigureWidget(fig)

    @render_widget
    def summary_drug_hist():
        if not raw_drugs or df is None: return go.FigureWidget()
        cpd_prom, prom_type_list = [], []
        
        for d in raw_drugs:
            log_col = f'log2 {d} R'
            if log_col in df.columns:
                valid = df[log_col].dropna()
                if len(valid) > 0:
                    cpd_prom.append((valid > 1).sum() / len(valid) * 100)
                    prom_type_list.append(type_dict.get(d, 'Unknown'))
                    
        prom_df = pd.DataFrame({'Promiscuity': cpd_prom, 'Type': prom_type_list})
        fig = px.histogram(prom_df, x='Promiscuity', color='Type', nbins=100, barmode='overlay', opacity=0.7, 
                           labels={'Promiscuity': 'Compound Promiscuity (% Sites Hit at R > 2)', 'Frequency': 'Frequency'})
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=40, r=40, t=10, b=40), legend_title_text='')
        return go.FigureWidget(fig)

    def get_site_compounds_bar_plot(site_label, sig_only):
        if not site_label or site_label == "Loading..." or df is None or df.empty: return go.FigureWidget()
        
        matching_rows = df[df['Labels'] == site_label]
        if matching_rows.empty: return go.FigureWidget()
        row = matching_rows.iloc[0]
        
        drug_data = []
        p_thresh = -np.log10(0.05)
        for d in raw_drugs:
            r_col = f'log2 {d} R'
            p_col = f'-log10 {d} p'
            if r_col in df.columns:
                r_val = row[r_col]
                if pd.isna(r_val): continue
                
                if sig_only and p_col in df.columns:
                    if pd.isna(row[p_col]) or row[p_col] <= p_thresh:
                        continue
                        
                drug_data.append({'Drug': d, 'R': 2**r_val})
                
        if not drug_data:
            return go.FigureWidget(px.bar(title="No compounds meet the criteria."))
            
        bar_df = pd.DataFrame(drug_data).sort_values('R', ascending=False).head(30)
        fig = px.bar(bar_df, x='Drug', y='R', color='R', color_continuous_scale='Reds')
        fig.update_layout(xaxis_title="", yaxis_title="R", plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=20, r=20, t=20, b=40), coloraxis_showscale=False)
        
        widget = go.FigureWidget(fig)
        if widget.data: 
            widget.data[0].on_click(lambda t, p, s: active_drug.set(t.x[p.point_inds[0]]) if p.point_inds else None)
        return widget

    @render_widget
    def summary_cancer_bar():
        return get_site_compounds_bar_plot(input.summary_cancer_site(), input.summary_sig_only())

    @render_widget
    def summary_ppi_bar():
        return get_site_compounds_bar_plot(input.summary_ppi_site(), input.summary_sig_only())

    # --- MOLECULE HTML RENDERERS ---
    def generate_molecule_html(drug):
        if not drug: return ui.p("Click a bar in the plots above to view the compound structure.", style="color: gray; font-style: italic;")
        smiles = smiles_dict.get(drug, "") 
        if not smiles: return ui.p(f"No SMILES available for {drug}.", style="color: gray; font-style: italic;")
        img_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(smiles)}/PNG"
        return ui.HTML(f'<div style="text-align: center;"><b>{drug}</b><br><img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px; margin-top: 5px;"></div>')

    @render.ui
    def molecule_ui_summary(): return generate_molecule_html(active_drug.get())

    @render.ui
    def molecule_ui_compound(): return generate_molecule_html(active_drug.get())

    @render.ui
    def molecule_ui_site(): return generate_molecule_html(active_drug.get())


    # --- COMPOUND & SITE VIEW PLOTS ---
    @reactive.Calc
    def available_pdbs():
        gene = input.target_gene()
        site_str = input.target_site_pos()
        if not gene or not site_str or gene == "No Data Found" or df is None or df.empty: return {}

        target = f"{gene}_Y{site_str}"
        matching_rows = df[df['Labels'] == target]
        if matching_rows.empty: return {}

        uniprot = str(matching_rows.iloc[0]['Protein Id']).split('|')[1].split('-')[0] if '|' in str(matching_rows.iloc[0]['Protein Id']) else str(matching_rows.iloc[0]['Protein Id']).split('-')[0]
        user_seq = str(matching_rows.iloc[0]['Sequence'])

        try:
            data = json.loads(urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{uniprot}.json").read().decode('utf-8'))
            pdb_list = []
            site_num = int(site_str) if str(site_str).isdigit() else -1

            for ref in data.get('uniProtKBCrossReferences', []):
                if ref['database'] == 'PDB':
                    pdb_id, method, res_val, has_site, coverage = ref['id'], "Unknown", 999.0, False, 0
                    for prop in ref.get('properties', []):
                        if prop['key'] == 'Method': method = prop['value']
                        elif prop['key'] == 'Resolution' and prop['value'] != '-':
                            try: res_val = float(prop['value'].replace('A', '').strip())
                            except: pass
                        elif prop['key'] == 'Chains' and '=' in prop['value']:
                            for start_str, end_str in re.findall(r'(\d+)-(\d+)', prop['value'].split('=', 1)[1]):
                                start, end = int(start_str), int(end_str)
                                coverage += (end - start + 1)
                                if start <= site_num <= end: has_site = True
                                
                    if has_site:
                        chain_id, mapped_site_pos = verify_and_map_site(pdb_id, site_str, user_seq)
                        if chain_id is None:
                            has_site = False

                    method_score = 1 if 'EM' in method.upper() else 2 if 'X-RAY' in method.upper() else 3 if 'NMR' in method.upper() else 4
                    pdb_list.append({'id': pdb_id, 'has_site': has_site, 'method_score': method_score, 'res_val': res_val, 'coverage': coverage, 'method_label': method, 'res_label': "N/A" if res_val == 999.0 else f"{res_val}Å"})
            
            pdb_list.sort(key=lambda x: (not x['has_site'], x['method_score'], x['res_val'], -x['coverage']))
            return {p['id']: f"{p['id']}|{'Site Present' if p['has_site'] else 'Site Absent'}|{p['method_label']}|{p['res_label']}|{p['coverage']}aa" for p in pdb_list}
        except Exception: return {}

    @reactive.Effect
    def update_pdb_dropdown():
        pdb_dict = available_pdbs()
        if pdb_dict: ui.update_select("pdb_selector", choices=pdb_dict, selected=list(pdb_dict.keys())[0])
        else: ui.update_select("pdb_selector", choices={"none": "No PDBs Found"}, selected="none")

    @reactive.Effect
    def update_ppi_dropdown():
        gene = input.target_gene()
        site_str = input.target_site_pos()
        if not gene or not site_str or ppi_df.empty:
            ui.update_select("ppi_selector", choices={"none": "No PPI Data Available"}, selected="none")
            return

        valid_ppis = ppi_df[ppi_df['Target'] == f"{gene}_Y{site_str}"].copy()
        if valid_ppis.empty: ui.update_select("ppi_selector", choices={"none": "No PPI interfaces found"}, selected="none")
        else:
            valid_ppis = valid_ppis.sort_values('Min_Distance', ascending=True)
            choices = {row['PDB_ID']: f"{row['PDB_ID']}|{row['Min_Distance']:.1f} Å" for _, row in valid_ppis.iterrows()}
            ui.update_select("ppi_selector", choices=choices, selected=list(choices.keys())[0])


    @reactive.Calc
    def plot_data():
        data_type, threshold, n_labels, color_mode = active_drug.get(), input.threshold(), input.n_labels(), input.color_mode()
        
        # Guard if no drug is selected yet (e.g. fresh launch on Compound tab)
        if not data_type or f'log2 {data_type} R' not in df.columns:
            return pd.DataFrame(), []

        plot_df = pd.DataFrame({
            'Site': np.arange(len(df[f'log2 {data_type} R'])), 'R_plot': 2**df[f'log2 {data_type} R'], 'R': 2**df[f'log2 {data_type} R'],
            'P-value': 10**-df[f'-log10 {data_type} p'], 'log2 R': df[f'log2 {data_type} R'], '-log10 p': df[f'-log10 {data_type} p'],
            'Label': df['Labels'], 'Gene Symbol': df['Gene Symbol'], 'ID': df['Protein Id'], 'Description': df['Description'],
        }).replace([np.inf, -np.inf], np.nan).dropna(subset=['R', 'P-value', 'log2 R', '-log10 p'])
        
        plot_df['Site Rank'] = plot_df['R'].rank(ascending=False).astype(int).astype(str)+'/'+str(len(plot_df))
        plot_df['% Site Rank'] = plot_df['R'].rank(ascending=False)/len(plot_df)*100
        plot_df['color'] = 'non-significant'

        if color_mode == "Above Threshold": plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
        elif color_mode == "Highlight Custom List": plot_df.loc[plot_df['Gene Symbol'].str.upper().isin([g.strip().upper() for g in input.custom_list().split(',') if g.strip()]), 'color'] = 'highlight'
        elif color_mode == "Cancer-Driver List": plot_df.loc[plot_df['Gene Symbol'].isin(cancer['Gene']), 'color'] = 'highlight'
        elif color_mode == "Sites with PPIs": plot_df.loc[plot_df['Label'].isin(set(ppi_df['Target'].dropna().unique())), 'color'] = 'highlight'
        else:
            plot_df.loc[(plot_df['R'] > threshold), 'color'] = 'high'
            plot_df['alpha'] = np.where(plot_df['R'] > threshold, 1.0, 0.2)

        sig_genes = plot_df[plot_df['color'] != 'non-significant'] if color_mode != "P-Value Gradient" else plot_df[plot_df['color'] == 'high']
        return plot_df, [idx for idx in sig_genes.nlargest(min(n_labels, len(sig_genes)), 'R').index if idx in plot_df.index]

    hover_dict = {'Site': False, 'color': False, 'R_plot': False, 'ID': True, 'Label': True, 'Description': True, 'log2 R': False, '-log10 p': False, 'R': ':.3f', 'P-value': ':.4f', 'Site Rank': True, '% Site Rank': ':.5f'}

    @render_widget
    def site_plot():
        plot_df, valid_indices = plot_data()
        if plot_df.empty: return go.FigureWidget()

        fig = px.scatter(plot_df, x='Site', y='R_plot', color='color', hover_data=hover_dict, color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD', 'highlight': '#D62728'}) if input.color_mode() != "P-Value Gradient" else px.scatter(plot_df, x='Site', y='R_plot', color='P-value', color_continuous_scale='Viridis_r', hover_data=hover_dict)
        if input.color_mode() == "P-Value Gradient": fig.update_traces(marker=dict(opacity=plot_df['alpha']))
        else: fig.update_traces(marker=dict(opacity=0.2 if input.color_mode() == "Above Threshold" else 0.15), selector=dict(name='non-significant'))
        
        for idx in valid_indices:
            r = plot_df.loc[idx]
            fig.add_annotation(x=r['Site'], y=r['R'], text=r['Label'], showarrow=True, arrowsize=1, arrowwidth=1, arrowcolor="#444", ax=20, ay=-30, font=dict(size=10, color="black"), bgcolor="white", bordercolor="#c7c7c7", borderwidth=1, borderpad=3)
        fig.add_hline(y=input.threshold(), line_width=1, line_dash='dash')
        fig.update_layout(yaxis_title="R", plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))
        widget = go.FigureWidget(fig)
        widget._config = {**(widget._config or {}), 'edits': {'annotationTail': True}}
        return widget

    @render_widget
    def volcano_plot():
        volcano_df, _ = plot_data() 
        if volcano_df.empty: return go.FigureWidget()

        safe_thresh = input.threshold() if input.threshold() > 0 else 1e-9 
        sig_mask = (volcano_df['log2 R'] > np.log2(safe_thresh)) & (volcano_df['-log10 p'] > -np.log10(0.05))
        volcano_df['color'] = 'non-significant'

        if input.color_mode() == "Highlight Custom List": volcano_df.loc[sig_mask & volcano_df['Gene Symbol'].str.upper().isin([g.strip().upper() for g in input.custom_list().split(',') if g.strip()]), 'color'] = 'highlight'
        elif input.color_mode() == "Cancer-Driver List": volcano_df.loc[sig_mask & volcano_df['Gene Symbol'].isin(cancer['Gene']), 'color'] = 'highlight'
        elif input.color_mode() == "Sites with PPIs": volcano_df.loc[sig_mask & volcano_df['Label'].isin(set(ppi_df['Target'].dropna().unique())), 'color'] = 'highlight'
        else: volcano_df.loc[sig_mask, 'color'] = 'high'

        colored_genes = volcano_df[volcano_df['color'] != 'non-significant']
        valid_indices = colored_genes.nlargest(min(input.n_labels(), len(colored_genes)), 'R').index.tolist() if not colored_genes.empty else []

        if input.color_mode() == "P-Value Gradient":
            fig = px.scatter(volcano_df[volcano_df['color'] == 'non-significant'], x='log2 R', y='-log10 p', color_discrete_sequence=['#dddddd'], hover_data=hover_dict)
            fig.update_traces(marker=dict(opacity=0.2), hovertemplate=None)
            sig_df = volcano_df[volcano_df['color'] == 'high']
            if not sig_df.empty:
                fig2 = px.scatter(sig_df, x='log2 R', y='-log10 p', color='P-value', color_continuous_scale='Viridis_r', hover_data=hover_dict)
                for t in fig2.data: fig.add_trace(t)
                fig.layout.coloraxis = fig2.layout.coloraxis
        else:
            fig = px.scatter(volcano_df, x='log2 R', y='-log10 p', color='color', color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD', 'highlight': '#D62728'}, hover_data=hover_dict)
            fig.update_traces(marker=dict(opacity=0.2), selector=dict(name='non-significant'))

        for idx in valid_indices:
            r = volcano_df.loc[idx]
            fig.add_annotation(x=r['log2 R'], y=r['-log10 p'], text=r['Label'], showarrow=True, arrowsize=1, arrowwidth=1, arrowcolor="#444", ax=20, ay=-30, font=dict(size=10, color="black"), bgcolor="white", bordercolor="#c7c7c7", borderwidth=1, borderpad=3)
        fig.add_vline(x=np.log2(safe_thresh), line_width=1, line_dash='dash')
        fig.add_hline(y=-np.log10(0.05), line_width=1, line_dash='dash')
        fig.update_layout(xaxis_title="log<sub>2</sub> Fold Change", yaxis_title="-log<sub>10</sub> P-Value", plot_bgcolor='white', paper_bgcolor='white', showlegend=False, margin=dict(l=40, r=40, t=10, b=40))
        widget = go.FigureWidget(fig)
        widget._config = {**(widget._config or {}), 'edits': {'annotationTail': True}}
        return widget

    @render_widget
    def site_profile_plot():
        gene, site = input.target_gene(), input.target_site_pos()
        if not gene or not site or df.empty: return go.FigureWidget(px.bar(title="No Site Selected"))
        matching_rows = df[df['Labels'] == f"{gene}_Y{site}"]
        if matching_rows.empty: return go.FigureWidget(px.bar(title="Loading..."))
             
        row = matching_rows.iloc[0]
        drug_data = [{'Drug': d, 'R': 2**row[f'log2 {d} R']} for d in raw_drugs if f'log2 {d} R' in df.columns]
        bar_df = pd.DataFrame(drug_data).dropna().sort_values('R', ascending=False)
        fig = px.bar(bar_df, x='Drug', y='R', color='R', color_continuous_scale='Reds', hover_data={'Drug': True, 'R': ':.3f'})
        fig.add_hline(y=input.threshold(), line_width=1, line_dash='dash', line_color='red')
        fig.update_layout(plot_bgcolor='white', paper_bgcolor='white', margin=dict(l=20, r=20, t=40, b=20), coloraxis_showscale=False)
        widget = go.FigureWidget(fig)
        if widget.data: widget.data[0].on_click(lambda t, p, s: active_drug.set(t.x[p.point_inds[0]]) if p.point_inds else None)
        return widget

    @render.ui
    def alphafold_viewer():
        gene, site_pos = input.target_gene(), input.target_site_pos()
        if not gene or not site_pos or df.empty: return ui.p("Please select a Gene and Site.")
        matching_rows = df[df['Labels'] == f"{gene}_Y{site_pos}"]
        if matching_rows.empty: return ui.p("Loading structure...")

        row = matching_rows.iloc[0]
        uniprot = str(row['Protein Id']).split('|')[1].split('-')[0] if '|' in str(row['Protein Id']) else str(row['Protein Id']).split('-')[0]

        chain_map_js = json.dumps({'A': str(row['Gene Symbol'])})
        hover_js = f"""function(atom,viewer,event,container) {{ if(!atom.label) {{ var cName = {chain_map_js}[atom.chain] || "Unknown Molecule"; atom.label = viewer.addLabel(atom.resn + " " + atom.resi + " (" + cName + ")", {{position: atom, backgroundColor: '#2b2b2b', fontColor: 'white', backgroundOpacity: 0.85, fontSize: 12, borderRadius: 5}}); }} }}"""

        try: pdb_data = urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v6.cif").read().decode('utf-8').replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception: return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve structure for {uniprot}.</div>")

        view = py3Dmol.view(width="100%", height=500)
        view.addModel(pdb_data, "cif")
        
        # Color by pLDDT (AlphaFold B-factor bins)
        view.setStyle({'b': [-100, 50]}, {'cartoon': {'color': '#ff7d45'}})   # Orange: Very low (<50)
        view.addStyle({'b': [50, 70]}, {'cartoon': {'color': '#ffe500'}})      # Yellow: Low (50-70)
        view.addStyle({'b': [70, 90]}, {'cartoon': {'color': '#65cbff'}})      # Cyan: Confident (70-90)
        view.addStyle({'b': [90, 1000]}, {'cartoon': {'color': '#0053d6'}})    # Blue: Very high (>90)
        
        # Highlight the specific site
        view.addStyle({'resi': site_pos, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
        
        view.zoomTo({'resi': site_pos})
        view.setHoverable({}, True, hover_js, unhover_js)

        return ui.HTML(f'''
        <div style="margin-bottom: 5px; font-size: 0.9em; line-height: 1.3;">
            <br><b>AlphaFold: <a href="https://alphafold.ebi.ac.uk/entry/AF-{uniprot}-F1" target="_blank">{uniprot}</a></b>
            <div style="margin-top: 5px; display: flex; gap: 10px; font-size: 0.85em; color: #555;">
                <span style="display: flex; align-items: center; gap: 4px;"><div style="width: 12px; height: 12px; background: #0053d6; border-radius: 2px;"></div> &gt;90 (Very High)</span>
                <span style="display: flex; align-items: center; gap: 4px;"><div style="width: 12px; height: 12px; background: #65cbff; border-radius: 2px;"></div> 70-90 (Confident)</span>
                <span style="display: flex; align-items: center; gap: 4px;"><div style="width: 12px; height: 12px; background: #ffe500; border-radius: 2px;"></div> 50-70 (Low)</span>
                <span style="display: flex; align-items: center; gap: 4px;"><div style="width: 12px; height: 12px; background: #ff7d45; border-radius: 2px;"></div> &lt;50 (Very Low)</span>
            </div>
        </div>
        <iframe src="data:text/html;base64,{base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')}" style="width: 100%; height: 500px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>
        ''')
    
    @render.ui
    def pdb_viewer():
        pdb_id, gene, site_pos = input.pdb_selector(), input.target_gene(), input.target_site_pos()
        if not pdb_id or pdb_id in ["Loading...", "none"]: return ui.HTML("<div style='color:orange'><b>Note:</b> No PDB structure found.</div>")
        if pdb_id not in available_pdbs(): return ui.p("Syncing...")

        matching_rows = df[df['Labels'] == f"{gene}_Y{site_pos}"]
        if matching_rows.empty: return ui.p("")
        
        row = matching_rows.iloc[0]
        chain_id, mapped_site_pos = verify_and_map_site(pdb_id, site_pos, str(row['Sequence']))

        mapping_notice = ""
        if chain_id is None:
            # FULL STRUCTURE FALLBACK
            mapping_notice = f"<br><span style='color: #d62728; font-size: 0.85em;'><i>Warning: Site not observed in this PDB, showing full structure.</i></span>"
        elif mapped_site_pos != site_pos:
            mapping_notice = f"<br><span style='color: #d62728; font-size: 0.85em;'><i>Sequence alignment matched site to author position {mapped_site_pos}.</i></span>"

        chain_map_js = json.dumps(get_chain_names(pdb_id))
        hover_js = f"""function(atom,viewer,event,container) {{ if(!atom.label) {{ var cName = {chain_map_js}[atom.chain] || "Unknown Molecule"; atom.label = viewer.addLabel(atom.resn + " " + atom.resi + " (" + cName + ")", {{position: atom, backgroundColor: '#2b2b2b', fontColor: 'white', backgroundOpacity: 0.85, fontSize: 12, borderRadius: 5}}); }} }}"""

        try: pdb_data = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.cif").read().decode('utf-8').replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception: return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve PDB {pdb_id}.</div>")

        view = py3Dmol.view(width="100%", height=480)
        view.addModel(pdb_data, "cif")
        view.setStyle({'cartoon': {'colorscheme': 'chain'}})
        
        if chain_id is not None:
            target_sel = {'resi': mapped_site_pos, 'chain': chain_id}
            view.addStyle({**{'not': {'atom': ['N', 'C', 'O', 'OXT']}}, **target_sel}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
            view.zoomTo(target_sel)
        else:
            view.zoomTo() # Fallback zoom to whole structure
            
        view.setHoverable({}, True, hover_js, unhover_js)
        
        return ui.HTML(f'''
        <div style="margin-bottom: 5px; font-size: 0.9em; line-height: 1.3;"><b>PDB ID: <a href="https://www.rcsb.org/structure/{pdb_id}" target="_blank">{pdb_id}</a></b>{mapping_notice}</div>
        <iframe src="data:text/html;base64,{base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')}" style="width: 100%; height: 440px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>
        ''')

    @render.ui
    def ppi_viewer():
        pdb_id, target_site = input.ppi_selector(), input.target_site_pos()
        if not pdb_id or pdb_id in ["none", "Loading..."]: return ui.div()

        match = ppi_df[(ppi_df['Target'] == f"{input.target_gene()}_Y{target_site}") & (ppi_df['PDB_ID'] == pdb_id)]
        viz_cutoff = float(match['Min_Distance'].iloc[0]) + 1.0 if not match.empty else 5.0

        matching_rows = df[df['Labels'] == f"{input.target_gene()}_Y{target_site}"]
        if matching_rows.empty: return ui.div()

        chain_id, mapped_site_pos = verify_and_map_site(pdb_id, target_site, str(matching_rows.iloc[0]['Sequence']))

        mapping_notice = ""
        if chain_id is None:
            # FULL STRUCTURE FALLBACK
            mapping_notice = f"<br><span style='color: #d62728; font-size: 0.85em;'><i>Warning: Site not observed in this PDB, showing full structure.</i></span>"
        elif mapped_site_pos != target_site:
            mapping_notice = f"<br><span style='color: #d62728; font-size: 0.85em;'><i>Sequence alignment matched site to author position {mapped_site_pos}.</i></span>"

        chain_map_js = json.dumps(get_chain_names(pdb_id))
        hover_js = f"""function(atom,viewer,event,container) {{ if(!atom.label) {{ var cName = {chain_map_js}[atom.chain] || "Unknown Molecule"; atom.label = viewer.addLabel(atom.resn + " " + atom.resi + " (" + cName + ")", {{position: atom, backgroundColor: '#2b2b2b', fontColor: 'white', backgroundOpacity: 0.85, fontSize: 12, borderRadius: 5}}); }} }}"""

        try: pdb_data = urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif").read().decode('utf-8').replace('\\', '\\\\').replace('`', '\\`').replace('$', '\\$')
        except Exception: return ui.HTML(f"<div style='color:red; padding:20px;'><b>Error:</b> Could not retrieve {pdb_id}.</div>")

        view = py3Dmol.view(width="100%", height=480)
        view.addModel(pdb_data, "cif")
        view.setStyle({'cartoon': {'colorscheme': 'chain'}})

        if chain_id is not None:
            target_sel = {'resi': mapped_site_pos, 'chain': chain_id}
            view.addStyle({'within': {'distance': viz_cutoff, 'sel': target_sel}, 'byres': True, 'not': {'atom': ['N', 'C', 'O', 'OXT']}}, {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.2}})
            view.addStyle({**{'not': {'atom': ['N', 'C', 'O', 'OXT']}}, **target_sel}, {'stick': {'colorscheme': 'redCarbon', 'radius': 0.2}})
            view.zoomTo(target_sel)
        else:
            view.zoomTo() # Fallback zoom to whole structure

        view.setHoverable({}, True, hover_js, unhover_js)

        return ui.HTML(f'''
        <div style="margin-bottom: 5px; font-size: 0.9em; line-height: 1.3;">
            <b>PDB ID: <a href="https://www.rcsb.org/structure/{pdb_id}" target="_blank">{pdb_id}</a></b><br>
            {mapping_notice}
        </div>
        <iframe src="data:text/html;base64,{base64.b64encode(view._make_html().encode('utf-8')).decode('utf-8')}" style="width: 100%; height: 440px; border: 1px solid #eee; border-radius: 5px; overflow: hidden;"></iframe>
        ''')

# --- 4. Run App ---
app = App(app_ui, server)