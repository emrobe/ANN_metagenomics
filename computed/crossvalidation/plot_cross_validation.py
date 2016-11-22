from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.plotly as py
import plotly.graph_objs as go

import pandas as pd

# Read data from a csv
z_data = pd.read_csv('surface_plot.csv')
print z_data.as_matrix()
data = [
    go.Surface(
        z=z_data.as_matrix(),
	x = [1,10,20,30,50,100],
	y = [1,5,10,20,40,60,100,150,200,400,600,800,1200,2000],
    )
]
layout = go.Layout(
    autosize=False,
#    scene=go.Scene(
#           font=dict(size=18),
#           xaxis=go.XAxis(title='Hidden nodes'),
#           yaxis=go.YAxis(title='Training epochs'),
#           zaxis=go.ZAxis(title='Accuracy')
#    ),
    scene = dict(
        xaxis = dict( title='Hidden nodes' ),
        yaxis = dict( title='Training epochs' ),
        zaxis = dict( title='Accuracy' ),
    ),
    xaxis=dict(title="Hidden Nodes"),
    yaxis=dict(title='Epochs Trained'),
    width=1600,
    height=1600,
    margin=dict(
        l=65,
        r=50,
        b=65,
        t=90
    ),
    font=dict(size=18),
)
fig = go.Figure(data=data, layout=layout)
plot(fig, filename='crossvalidation')
