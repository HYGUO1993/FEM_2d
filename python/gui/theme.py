"""
Theme for FEM_2d GUI Application

Provides dark theme CSS (QSS), color palettes, and matplotlib configuration.
"""

from PySide6.QtGui import QPalette, QColor

# Catppuccin Macchiato inspired palette
PALETTE = {
    "base": "#24273A",
    "mantle": "#1E2030",
    "crust": "#181926",
    "text": "#CAD3F5",
    "subtext0": "#A5ADCB",
    "subtext1": "#B8C0E0",
    "surface0": "#363A4F",
    "surface1": "#494D64",
    "surface2": "#5B6078",
    "overlay0": "#6E738D",
    "overlay1": "#8087A2",
    "overlay2": "#939AB7",
    "blue": "#8AADF4",
    "lavender": "#B7BDF8",
    "sapphire": "#7DC4E4",
    "sky": "#91D7E3",
    "teal": "#8BD5CA",
    "green": "#A6DA95",
    "yellow": "#EED49F",
    "peach": "#F5A97F",
    "maroon": "#EE99A0",
    "red": "#ED8796",
    "mauve": "#C6A0F6",
    "pink": "#F5BDE6",
    "flamingo": "#F0C6C6",
    "rosewater": "#F4DBD6",
}

DARK_QSS = f"""
QWidget {{
    background-color: {PALETTE['base']};
    color: {PALETTE['text']};
    font-family: "Segoe UI", "Microsoft YaHei", sans-serif;
    font-size: 10pt;
}}

QMainWindow {{
    background-color: {PALETTE['crust']};
}}

QMenuBar {{
    background-color: {PALETTE['crust']};
    color: {PALETTE['text']};
    border-bottom: 1px solid {PALETTE['surface0']};
}}

QMenuBar::item:selected {{
    background-color: {PALETTE['surface0']};
}}

QMenu {{
    background-color: {PALETTE['mantle']};
    border: 1px solid {PALETTE['surface0']};
}}

QMenu::item:selected {{
    background-color: {PALETTE['surface1']};
}}

QToolBar {{
    background-color: {PALETTE['mantle']};
    border-bottom: 1px solid {PALETTE['surface0']};
    spacing: 5px;
}}

QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: 4px;
    padding: 3px;
}}

QToolButton:hover {{
    background-color: {PALETTE['surface1']};
}}

QToolButton:pressed {{
    background-color: {PALETTE['surface2']};
}}

QDockWidget {{
    background-color: {PALETTE['mantle']};
    color: {PALETTE['blue']};
    font-weight: bold;
}}

QDockWidget::title {{
    background-color: {PALETTE['crust']};
    padding: 6px;
}}

QTreeView, QTreeWidget, QTextEdit {{
    background-color: {PALETTE['mantle']};
    border: 1px solid {PALETTE['surface0']};
    border-radius: 4px;
}}

QTreeView::item:selected, QTreeWidget::item:selected {{
    background-color: {PALETTE['surface1']};
    color: {PALETTE['blue']};
}}

QScrollBar:vertical {{
    background-color: {PALETTE['mantle']};
    width: 12px;
    margin: 0px;
}}

QScrollBar::handle:vertical {{
    background-color: {PALETTE['surface2']};
    min-height: 20px;
    border-radius: 6px;
    margin: 2px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {PALETTE['overlay0']};
}}

QScrollBar:horizontal {{
    background-color: {PALETTE['mantle']};
    height: 12px;
    margin: 0px;
}}

QScrollBar::handle:horizontal {{
    background-color: {PALETTE['surface2']};
    min-width: 20px;
    border-radius: 6px;
    margin: 2px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {PALETTE['overlay0']};
}}

QPushButton {{
    background-color: {PALETTE['surface1']};
    color: {PALETTE['text']};
    border: 1px solid {PALETTE['surface2']};
    border-radius: 4px;
    padding: 5px 15px;
}}

QPushButton:hover {{
    background-color: {PALETTE['surface2']};
    border: 1px solid {PALETTE['overlay0']};
}}

QPushButton:pressed {{
    background-color: {PALETTE['surface0']};
}}

QDoubleSpinBox, QSpinBox {{
    background-color: {PALETTE['crust']};
    border: 1px solid {PALETTE['surface1']};
    border-radius: 4px;
    padding: 4px;
}}

QDoubleSpinBox:hover, QSpinBox:hover {{
    border: 1px solid {PALETTE['blue']};
}}

QCheckBox {{
    spacing: 8px;
}}

QCheckBox::indicator {{
    width: 18px;
    height: 18px;
    border-radius: 4px;
    border: 1px solid {PALETTE['surface2']};
    background-color: {PALETTE['crust']};
}}

QCheckBox::indicator:hover {{
    border: 1px solid {PALETTE['blue']};
}}

QCheckBox::indicator:checked {{
    background-color: {PALETTE['blue']};
    border: 1px solid {PALETTE['blue']};
    image: url(check.png); /* A small white checkmark would go here if provided, otherwise relies on OS styling usually or requires base64 */
}}

QStatusBar {{
    background-color: {PALETTE['crust']};
    border-top: 1px solid {PALETTE['surface0']};
}}
"""

def apply_theme(app):
    """Apply the dark theme to the QApplication instance."""
    app.setStyleSheet(DARK_QSS)
    
    # Also set dark palette as fallback
    palette = QPalette()
    palette.setColor(QPalette.Window, QColor(PALETTE['base']))
    palette.setColor(QPalette.WindowText, QColor(PALETTE['text']))
    palette.setColor(QPalette.Base, QColor(PALETTE['mantle']))
    palette.setColor(QPalette.AlternateBase, QColor(PALETTE['surface0']))
    palette.setColor(QPalette.ToolTipBase, QColor(PALETTE['crust']))
    palette.setColor(QPalette.ToolTipText, QColor(PALETTE['text']))
    palette.setColor(QPalette.Text, QColor(PALETTE['text']))
    palette.setColor(QPalette.Button, QColor(PALETTE['surface1']))
    palette.setColor(QPalette.ButtonText, QColor(PALETTE['text']))
    palette.setColor(QPalette.BrightText, QColor(PALETTE['red']))
    palette.setColor(QPalette.Link, QColor(PALETTE['blue']))
    palette.setColor(QPalette.Highlight, QColor(PALETTE['surface2']))
    palette.setColor(QPalette.HighlightedText, QColor(PALETTE['text']))
    
    app.setPalette(palette)


def setup_matplotlib_dark_theme():
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')
    
    import matplotlib as mpl
    mpl.rcParams.update({
        "figure.facecolor": PALETTE['crust'],
        "axes.facecolor": PALETTE['base'],
        "axes.edgecolor": PALETTE['surface2'],
        "axes.labelcolor": PALETTE['text'],
        "text.color": PALETTE['text'],
        "xtick.color": PALETTE['subtext0'],
        "ytick.color": PALETTE['subtext0'],
        "grid.color": PALETTE['surface0'],
        "grid.alpha": 0.5,
        "grid.linestyle": "--",
        "legend.facecolor": PALETTE['mantle'],
        "legend.edgecolor": PALETTE['surface1'],
        "lines.color": PALETTE['text'],
        "patch.edgecolor": PALETTE['surface1'],
        "savefig.facecolor": PALETTE['crust'],
        "savefig.edgecolor": PALETTE['crust']
    })

