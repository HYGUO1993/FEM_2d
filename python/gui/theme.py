"""
Theme for FEM_2d GUI Application

Provides dark theme CSS (QSS), color palettes, and matplotlib configuration.
"""

from PyQt6.QtGui import QPalette, QColor

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

# Accent colour used for primary action buttons and active states
ACCENT = PALETTE["blue"]

DARK_QSS = f"""
/* ── Global ───────────────────────────────────────────────────── */
QWidget {{
    background-color: {PALETTE['base']};
    color: {PALETTE['text']};
    font-family: "Segoe UI", "Microsoft YaHei", "PingFang SC", sans-serif;
    font-size: 10pt;
}}

QMainWindow {{
    background-color: {PALETTE['crust']};
}}

/* ── Menu bar ─────────────────────────────────────────────────── */
QMenuBar {{
    background-color: {PALETTE['crust']};
    color: {PALETTE['text']};
    border-bottom: 1px solid {PALETTE['surface0']};
    padding: 2px 4px;
}}

QMenuBar::item {{
    padding: 4px 10px;
    border-radius: 4px;
}}

QMenuBar::item:selected {{
    background-color: {PALETTE['surface1']};
}}

QMenu {{
    background-color: {PALETTE['mantle']};
    border: 1px solid {PALETTE['surface1']};
    border-radius: 6px;
    padding: 4px 0;
}}

QMenu::item {{
    padding: 6px 24px 6px 16px;
    border-radius: 4px;
    margin: 1px 4px;
}}

QMenu::item:selected {{
    background-color: {PALETTE['surface1']};
    color: {PALETTE['blue']};
}}

QMenu::separator {{
    height: 1px;
    background-color: {PALETTE['surface0']};
    margin: 4px 8px;
}}

/* ── Toolbar ──────────────────────────────────────────────────── */
QToolBar {{
    background-color: {PALETTE['mantle']};
    border-bottom: 1px solid {PALETTE['surface0']};
    spacing: 4px;
    padding: 3px 6px;
}}

QToolBar::separator {{
    width: 1px;
    background-color: {PALETTE['surface1']};
    margin: 4px 6px;
}}

QToolButton {{
    background-color: transparent;
    border: 1px solid transparent;
    border-radius: 6px;
    padding: 5px 8px;
    color: {PALETTE['text']};
}}

QToolButton:hover {{
    background-color: {PALETTE['surface1']};
    border-color: {PALETTE['surface2']};
}}

QToolButton:pressed {{
    background-color: {PALETTE['surface0']};
}}

QToolButton:checked {{
    background-color: {PALETTE['surface1']};
    border-color: {PALETTE['blue']};
    color: {PALETTE['blue']};
}}

/* ── Side-panel tab bar ───────────────────────────────────────── */
QTabBar {{
    background: transparent;
}}

QTabBar::tab {{
    background-color: transparent;
    color: {PALETTE['subtext0']};
    border: none;
    padding: 8px 16px;
    font-size: 9pt;
    border-bottom: 2px solid transparent;
}}

QTabBar::tab:selected {{
    color: {PALETTE['blue']};
    border-bottom: 2px solid {PALETTE['blue']};
    font-weight: bold;
}}

QTabBar::tab:hover:!selected {{
    color: {PALETTE['text']};
    background-color: {PALETTE['surface0']};
    border-radius: 4px 4px 0 0;
}}

QTabWidget::pane {{
    border: none;
    background-color: {PALETTE['mantle']};
}}

/* ── Tree / List ──────────────────────────────────────────────── */
QTreeView, QTreeWidget, QListWidget {{
    background-color: {PALETTE['mantle']};
    border: 1px solid {PALETTE['surface0']};
    border-radius: 6px;
    alternate-background-color: {PALETTE['base']};
    outline: none;
}}

QTreeView::item, QTreeWidget::item, QListWidget::item {{
    padding: 3px 6px;
    border-radius: 4px;
}}

QTreeView::item:selected, QTreeWidget::item:selected, QListWidget::item:selected {{
    background-color: {PALETTE['surface1']};
    color: {PALETTE['blue']};
}}

QTreeView::item:hover, QTreeWidget::item:hover, QListWidget::item:hover {{
    background-color: {PALETTE['surface0']};
}}

QHeaderView::section {{
    background-color: {PALETTE['crust']};
    color: {PALETTE['subtext0']};
    border: none;
    padding: 4px 8px;
    font-weight: bold;
    font-size: 9pt;
}}

/* ── Scroll bars ──────────────────────────────────────────────── */
QScrollBar:vertical {{
    background-color: transparent;
    width: 10px;
    margin: 0;
}}

QScrollBar::handle:vertical {{
    background-color: {PALETTE['surface2']};
    min-height: 24px;
    border-radius: 5px;
    margin: 2px 2px;
}}

QScrollBar::handle:vertical:hover {{
    background-color: {PALETTE['overlay1']};
}}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}

QScrollBar:horizontal {{
    background-color: transparent;
    height: 10px;
    margin: 0;
}}

QScrollBar::handle:horizontal {{
    background-color: {PALETTE['surface2']};
    min-width: 24px;
    border-radius: 5px;
    margin: 2px 2px;
}}

QScrollBar::handle:horizontal:hover {{
    background-color: {PALETTE['overlay1']};
}}

QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}

/* ── Text editor ──────────────────────────────────────────────── */
QPlainTextEdit, QTextEdit {{
    background-color: {PALETTE['mantle']};
    border: 1px solid {PALETTE['surface0']};
    border-radius: 6px;
    selection-background-color: {PALETTE['surface2']};
    selection-color: {PALETTE['text']};
}}

QPlainTextEdit:focus, QTextEdit:focus {{
    border-color: {PALETTE['blue']};
}}

/* ── Buttons ──────────────────────────────────────────────────── */
QPushButton {{
    background-color: {PALETTE['surface0']};
    color: {PALETTE['text']};
    border: 1px solid {PALETTE['surface1']};
    border-radius: 6px;
    padding: 6px 18px;
    font-weight: 500;
}}

QPushButton:hover {{
    background-color: {PALETTE['surface1']};
    border-color: {PALETTE['overlay0']};
}}

QPushButton:pressed {{
    background-color: {PALETTE['surface0']};
    border-color: {PALETTE['blue']};
}}

QPushButton:disabled {{
    color: {PALETTE['overlay0']};
    background-color: {PALETTE['mantle']};
    border-color: {PALETTE['surface0']};
}}

/* Primary accent button */
QPushButton[accent="true"] {{
    background-color: {PALETTE['blue']};
    color: {PALETTE['crust']};
    border: none;
    font-weight: bold;
}}

QPushButton[accent="true"]:hover {{
    background-color: {PALETTE['lavender']};
}}

QPushButton[accent="true"]:pressed {{
    background-color: {PALETTE['sapphire']};
}}

/* ── Spin boxes ───────────────────────────────────────────────── */
QDoubleSpinBox, QSpinBox {{
    background-color: {PALETTE['crust']};
    border: 1px solid {PALETTE['surface1']};
    border-radius: 5px;
    padding: 4px 6px;
    selection-background-color: {PALETTE['surface2']};
}}

QDoubleSpinBox:focus, QSpinBox:focus {{
    border-color: {PALETTE['blue']};
}}

QDoubleSpinBox::up-button, QSpinBox::up-button,
QDoubleSpinBox::down-button, QSpinBox::down-button {{
    background-color: {PALETTE['surface0']};
    border: none;
    border-radius: 3px;
    width: 16px;
}}

QDoubleSpinBox::up-button:hover, QSpinBox::up-button:hover,
QDoubleSpinBox::down-button:hover, QSpinBox::down-button:hover {{
    background-color: {PALETTE['surface2']};
}}

/* ── Checkboxes ───────────────────────────────────────────────── */
QCheckBox {{
    spacing: 8px;
    color: {PALETTE['text']};
}}

QCheckBox::indicator {{
    width: 16px;
    height: 16px;
    border-radius: 4px;
    border: 1.5px solid {PALETTE['surface2']};
    background-color: {PALETTE['crust']};
}}

QCheckBox::indicator:hover {{
    border-color: {PALETTE['blue']};
}}

QCheckBox::indicator:checked {{
    background-color: {PALETTE['blue']};
    border-color: {PALETTE['blue']};
}}

/* ── Labels ───────────────────────────────────────────────────── */
QLabel {{
    background: transparent;
}}

QLabel[heading="true"] {{
    color: {PALETTE['blue']};
    font-weight: bold;
    font-size: 10pt;
}}

QLabel[subheading="true"] {{
    color: {PALETTE['subtext0']};
    font-size: 9pt;
}}

/* ── Splitter ─────────────────────────────────────────────────── */
QSplitter::handle {{
    background-color: {PALETTE['surface0']};
}}

QSplitter::handle:horizontal {{
    width: 2px;
}}

QSplitter::handle:vertical {{
    height: 2px;
}}

QSplitter::handle:hover {{
    background-color: {PALETTE['blue']};
}}

/* ── Status bar ───────────────────────────────────────────────── */
QStatusBar {{
    background-color: {PALETTE['crust']};
    border-top: 1px solid {PALETTE['surface0']};
    font-size: 9pt;
    color: {PALETTE['subtext0']};
}}

QStatusBar::item {{
    border: none;
}}

/* ── Tooltip ──────────────────────────────────────────────────── */
QToolTip {{
    background-color: {PALETTE['mantle']};
    color: {PALETTE['text']};
    border: 1px solid {PALETTE['surface1']};
    border-radius: 4px;
    padding: 4px 8px;
}}

/* ── Frame / group box ────────────────────────────────────────── */
QGroupBox {{
    border: 1px solid {PALETTE['surface1']};
    border-radius: 6px;
    margin-top: 1em;
    padding-top: 8px;
    font-weight: bold;
    color: {PALETTE['subtext1']};
}}

QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 6px;
    left: 10px;
}}

/* ── Form layout label column ─────────────────────────────────── */
QFormLayout QLabel {{
    color: {PALETTE['subtext0']};
    font-size: 9pt;
}}
"""


def apply_theme(app):
    """Apply the dark theme to the QApplication instance."""
    app.setStyleSheet(DARK_QSS)

    palette = QPalette()
    palette.setColor(QPalette.ColorRole.Window, QColor(PALETTE['base']))
    palette.setColor(QPalette.ColorRole.WindowText, QColor(PALETTE['text']))
    palette.setColor(QPalette.ColorRole.Base, QColor(PALETTE['mantle']))
    palette.setColor(QPalette.ColorRole.AlternateBase, QColor(PALETTE['surface0']))
    palette.setColor(QPalette.ColorRole.ToolTipBase, QColor(PALETTE['crust']))
    palette.setColor(QPalette.ColorRole.ToolTipText, QColor(PALETTE['text']))
    palette.setColor(QPalette.ColorRole.Text, QColor(PALETTE['text']))
    palette.setColor(QPalette.ColorRole.Button, QColor(PALETTE['surface1']))
    palette.setColor(QPalette.ColorRole.ButtonText, QColor(PALETTE['text']))
    palette.setColor(QPalette.ColorRole.BrightText, QColor(PALETTE['red']))
    palette.setColor(QPalette.ColorRole.Link, QColor(PALETTE['blue']))
    palette.setColor(QPalette.ColorRole.Highlight, QColor(PALETTE['surface2']))
    palette.setColor(QPalette.ColorRole.HighlightedText, QColor(PALETTE['text']))

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
        "savefig.edgecolor": PALETTE['crust'],
    })

