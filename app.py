from shiny import App, ui, render
from shinywidgets import output_widget, render_widget
import pandas as pd
import numpy as np
import plotly.express as px

# --- 1. Mock Data Setup (Replace with your actual 'test' and 'drugs' data) ---
np.random.seed(42)
drugs = ['DrugA', 'DrugB', 'DrugC', 'DrugD', 'DrugE', 'DrugF', 'DrugG']
n_sites = 300

test = pd.DataFrame({
    'Labels': [f'GENE_{i}' for i in range(n_sites)],
    'Protein Id': [f'P{i:04d}' for i in range(n_sites)]
})

# Generate dummy columns for each drug to mimic your dataset
for drug in drugs:
    test[f'log2 {drug} R'] = np.random.normal(0.5, 1.2, n_sites)
    test[f'-log10 {drug} p'] = np.random.uniform(0, 3, n_sites)

# --- 2. Shiny UI ---
app_ui = ui.page_fluid(
    ui.h2("Site-Level Significance Plot"),
    ui.layout_sidebar(
        ui.sidebar(
            ui.input_select("data_type", "Select Drug:", choices=drugs, selected=drugs[6]),
            ui.input_numeric("threshold", "Log2 R Threshold:", value=2.0, step=0.5),
            ui.input_numeric("n_labels", "Number of Top Labels:", value=5, min=0, max=20),
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
    
    @render_widget
    def site_plot():
        # Get reactive inputs
        data_type = input.data_type()
        threshold = input.threshold()
        n_labels = input.n_labels()
        df = test 

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
        fig.add_hline(y=1, line_width=1, line_dash='dash')

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