import dash
import dash_bootstrap_components as dbc
import dash_mantine_components as dmc
from dash_extensions.enrich import html, DashProxy, LogTransform, ServersideOutputTransform, MultiplexerTransform
import feffery_antd_components as fac
from dash import DiskcacheManager
import diskcache
from dash_extensions.enrich import FileSystemBackend
from dash_iconify import DashIconify

background_callback_manager = DiskcacheManager(diskcache.Cache("/data1/share/omics-viewer/cache"))

dash._dash_renderer._set_react_version('18.2.0')
dbc_css = "/data1/share/omics-viewer/main/dbc.min.css"
app = DashProxy(
  __name__, 
  external_stylesheets=[
    dbc.themes.BOOTSTRAP,
    dmc.styles.DATES,
    dmc.styles.CODE_HIGHLIGHT,
    dmc.styles.CHARTS,
    dmc.styles.CAROUSEL,
    dmc.styles.NOTIFICATIONS,
    dmc.styles.NPROGRESS,
  ],
  external_scripts = [
    {'src': 'https://deno.land/x/corejs@v3.31.1/index.js', 'type': 'module'}
  ],
  transforms=[
    LogTransform(), 
    ServersideOutputTransform(backends=[FileSystemBackend("/data1/share/omics-viewer/cache")]), 
    MultiplexerTransform()
  ],
  use_pages=True,
  requests_pathname_prefix='/',
  background_callback_manager=background_callback_manager
)

header = dbc.Navbar(
  dbc.Container([
    dbc.Row([
      dbc.Col(dbc.NavbarBrand("SEU-3D", href="/", style={'color': 'white'})),
      dbc.Col(
        dmc.NavLink(
          label = "Documentation", 
          href="https://rainyblue-w.github.io/SEU-3D/",
          rightSection = DashIconify(icon='hugeicons:google-doc', width=24),
          variant='filled', autoContrast=True, color='dark', active=True,
        ),
      )
    ]),
    dbc.DropdownMenu(
        children=[
            dbc.DropdownMenuItem('spatial', href='/'),
            dbc.DropdownMenuItem("atlas", href='/atlasall'),
            dbc.DropdownMenuItem("atlas-cp", href='/atlascp'),
            dbc.DropdownMenuItem("atlas-endoderm", href='/atlased'),
            dbc.DropdownMenuItem("reik", href='/reik'),
            dbc.DropdownMenuItem("seqFISH", href='/seqfish'),
        ],
        nav=True,
        in_navbar=True,
        label="Dataset",
        style={'color': 'white'}
    ),
  ]),
  color="#2e2e2e",
  className='dbc-Navbar-main',
)

app.layout = dmc.MantineProvider(
  html.Div([
    header,
    dash.page_container
  ])
)


#run server
if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser()
  parser.add_argument('--port', type=int, default=8050)
  parser.add_argument('--debug', action='store_true')
  args = parser.parse_args()
  app.run(
    host='::',
    port=args.port,
    threaded=True,
    proxy=None,
    debug = args.debug,
    use_reloader = False
  )