import numpy as np
import pandas as pd
import urllib.parse
from pathlib import Path
import plotly.express as px
from shiny import App, ui, render
from shinywidgets import output_widget, render_widget
from functools import reduce

# 1. Setup Paths
app_dir = Path(__file__).parent
data_dir = app_dir / "data"

# 2. Safely load the combined CSVs
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
else:
    df = pd.DataFrame()
    drugs = ["No Data Found"]

# 3. Safely load the Smiles dictionary & Cancer List
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
    cancer = pd.DataFrame(columns=['Gene']) # Empty fallback

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Compound Engagement Plot"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_select("data_type", "Select Drug:", choices=drugs, selected=drugs[0] if drugs else None),
            ui.input_numeric("threshold", "R Threshold:", value=2.0, step=0.5),
            ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
            
            # --- COLORING OPTIONS ---
            ui.hr(),
            ui.input_select(
                "color_mode", "Color Points By:", 
                choices=["Above Threshold", "P-Value Gradient", "Cancer-Driver List", "Highlight Custom List"]
            ),
            ui.input_text("custom_list", "Genes to Highlight (comma-separated):", placeholder="e.g. MAPK1, EGFR"),
            
            # --- STRUCTURE UI ---
            ui.hr(),
            ui.h5("Chemical Structure"),
            ui.output_ui("molecule_ui"), 
            
            width=300
        ),
        ui.card(
            output_widget("site_plot"),
            full_screen=True
        )
    )
)

# --- 3. Shiny Server ---
def server(input, output, session):
    
    @render.ui
    def molecule_ui():
        data_type = input.data_type()
        smiles = smiles_dict.get(data_type, "") 
        
        if not smiles:
            return ui.p("No SMILES available for this drug.", style="color: gray; font-style: italic;")

        safe_smiles = urllib.parse.quote(smiles)
        img_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/image"
        return ui.HTML(f'<img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Chemical Structure">')

    @render_widget
    def site_plot():
        data_type = input.data_type()
        threshold = input.threshold()
        n_labels = input.n_labels()
        color_mode = input.color_mode()
        custom_list_text = input.custom_list()

        # Create DataFrame for plotting
        plot_df = pd.DataFrame({
            'Site': np.arange(len(df[f'log2 {data_type} R'])),
            'R': 2**df[f'log2 {data_type} R'],
            'p': 10**-df[f'-log10 {data_type} p'],
            'Label': df['Labels'],
            'Gene Symbol': df['Gene Symbol'],
            'ID': df['Protein Id'] + '|' + df['Labels'].str.split('_', expand=True).iloc[:, 1]
        })
        
        # Clean Data
        plot_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        plot_df.dropna(subset=['R', 'p'], inplace=True)

        # --- COLORING LOGIC ---
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


        # --- ANNOTATION (LABELING) LOGIC ---
        pos_count = min(n_labels, len(sig_genes))
        if pos_count > 0:
            top_positive = sig_genes.nlargest(pos_count, 'R').index
            valid_indices = [idx for idx in top_positive if idx in plot_df.index]
            
            for idx in valid_indices:
                row = plot_df.loc[idx]
                fig.add_annotation(
                    x=row['Site'], y=row['R'],
                    text=row['Label'],
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=1,
                    arrowwidth=1.5,
                    arrowcolor="#444",
                    ax=20,  # Shifts the label 20 pixels to the right
                    ay=-30, # Shifts the label 30 pixels up
                    font=dict(size=10, color="black"),
                    bgcolor="white",
                    bordercolor="#c7c7c7",
                    borderwidth=1,
                    borderpad=3
                )

        # --- FORMATTING ---
        fig.add_hline(y=threshold, line_width=1, line_dash='dash')

        fig.update_traces(marker=dict(size=6))

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            showlegend=False,
            margin=dict(l=40, r=40, t=40, b=40)
        )

        return fig

# --- 4. Run App ---
app = App(app_ui, server)