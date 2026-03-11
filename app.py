import numpy as np
import pandas as pd
import urllib.parse
from pathlib import Path
import plotly.express as px
from shiny import App, ui, render
from shinywidgets import output_widget, render_widget

# 1. Get the directory where app.py is currently running
app_dir = Path(__file__).parent

# --- 1. Data Setup ---
df = pd.read_csv(app_dir/"drug_R.csv")
drugs = list(set([i.split(' ')[1] for i in df.columns if 'log' in i]))

ypt_lib = pd.read_csv(app_dir/'ypt_library.csv')
smiles_dict = dict(zip(ypt_lib['Catalog ID'], ypt_lib['Smiles']))

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Site-Level Significance Plot"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_select("data_type", "Select Drug:", choices=drugs, selected=drugs[0]),
            ui.input_numeric("threshold", "Log2 R Threshold:", value=2.0, step=0.5),
            ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
            
            # --- NEW UI BOX FOR STRUCTURE ---
            ui.hr(),
            ui.h5("Chemical Structure"),
            ui.output_ui("molecule_ui"), # Outputs the HTML image
            
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
    
    # --- NEW SERVER FUNCTION FOR STRUCTURE ---
    @render.ui
    def molecule_ui():
        data_type = input.data_type()
        
        # Get the SMILES string
        smiles = smiles_dict.get(data_type, "") 
        
        if not smiles:
            return ui.p("No SMILES available for this drug.", style="color: gray; font-style: italic;")

        # URL-encode the SMILES string so special characters (+, #, /, etc.) don't break the web link
        safe_smiles = urllib.parse.quote(smiles)
        
        # Use the NCI CACTUS API to generate the image automatically
        img_url = f"https://cactus.nci.nih.gov/chemical/structure/{safe_smiles}/image"
        
        # Return the HTML image tag directly pointing to the API
        return ui.HTML(f'<img src="{img_url}" style="width:100%; max-width:250px; background-color:white; border: 1px solid #ddd; border-radius: 4px; padding: 5px;" alt="Chemical Structure">')

    @render_widget
    def site_plot():
        # Get reactive inputs
        data_type = input.data_type()
        threshold = input.threshold()
        n_labels = input.n_labels()

        # Create DataFrame for plotting
        plot_df = pd.DataFrame({
            'Site': np.arange(len(df[f'log2 {data_type} R'])),
            'log2 R': df[f'log2 {data_type} R'],
            '-log10 p': df[f'-log10 {data_type} p'],
            'Label': df['Labels'],
            'ID': df['Protein Id'] + '|' + df['Labels'].str.split('_', expand=True).iloc[:, 1]
        })

        # Define Significance Colors
        plot_df['color'] = 'non-significant'
        plot_df.loc[(plot_df['log2 R'] > np.log2(threshold)) & (plot_df['-log10 p'] > -np.log10(0.05)), 'color'] = 'high'

        # Labeling Logic
        plot_df['label'] = '' # Initialize with empty strings
        sig_genes = plot_df[plot_df['color'] == 'high']

        # Safe count handling for labels
        pos_count = min(n_labels, sum(plot_df['color'] == 'high'))
        
        if pos_count > 0:
            top_positive = sig_genes.nlargest(pos_count, 'log2 R').index
            indices_to_label = top_positive
            plot_df.loc[indices_to_label, 'label'] = plot_df.loc[indices_to_label, 'Label']

        # Plotting
        fig = px.scatter(
            plot_df,
            x='Site',
            y='log2 R',
            color='color',
            color_discrete_map={'non-significant': '#dddddd', 'high': '#4470AD', 'low': '#B95756'},
            text='label',  
            hover_data=['ID', 'log2 R', '-log10 p']
        )

        # Add thresholds
        fig.add_hline(y=threshold, line_width=1, line_dash='dash')

        # Formatting
        fig.update_traces(
            marker=dict(size=6),
            textposition='top center',
            textfont=dict(size=10, color='black')
        )

        fig.update_layout(
            plot_bgcolor='white',
            paper_bgcolor='white',
            showlegend=False,
            margin=dict(l=40, r=40, t=40, b=40)
        )

        return fig

# --- 4. Run App ---
app = App(app_ui, server)
