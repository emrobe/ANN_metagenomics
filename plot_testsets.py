from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot
import plotly.graph_objs as go
#from plotly.graph_objs import Bar, Scatter, Figure, Layout
import plotly.plotly as py
# Generate the figure

def testset_plot(x,y,z):

	annotations = []
	for n, row in enumerate(z):
	    for m, val in enumerate(row):
	        #var = z[n][m]
	        annotations.append(
	            dict(
	                text=str(round(val,3)),
	                x=x[m], y=y[n],
	                xref='x1', yref='y1',
	                font=dict(color='white' if val > 0.5 else 'black', size=20),
	                showarrow=False)
	            )

	colorscale = [[0, '#3D9970'], [1, '#001f3f']]  # custom colorscale
	trace = go.Heatmap(x=x, y=y, z=z, colorscale=colorscale, showscale=True)

	fig = go.Figure(data=[trace])
	fig['layout'].update(
	    annotations=annotations,
	    xaxis=dict(ticks=' ', ticksuffix='  ', side='top'),
	    # ticksuffix is a workaround to add a bit of padding
	    yaxis=dict(ticks=' ', ticksuffix='  '),
	    margin=go.Margin(l=200, t=200),
	    width=1400,
	    height=1400,
	    font=dict(size=20),
	    autosize=False
	)
	plot(fig, filename='result.html')
