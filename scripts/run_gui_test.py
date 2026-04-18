import sys
import os
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import QTimer

# Add project root and python dir to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))

from fem_gui import check_dependencies

def test_gui_screenshot():
    from gui.main_window import MainWindow
    from gui.theme import apply_theme
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    apply_theme(app)
    
    window = MainWindow()
    window.show()
    
    def simulate_user_actions():
        # Open default example
        window.open_default_example()
        
        # Trigger solution
        window.solve_model()
        
        # Give it a tiny bit of time to render
        QTimer.singleShot(1000, capture_and_exit)
        
    def capture_and_exit():
        # Grab screenshot of the main window
        pixmap = window.grab()
        screenshot_path = os.path.join(window.project_root, "build", "gui_screenshot.png")
        pixmap.save(screenshot_path)
        print(f"Saved screenshot to {screenshot_path}")
        app.quit()
        
    # Start simulating after the window is shown
    QTimer.singleShot(500, simulate_user_actions)
    
    sys.exit(app.exec())

if __name__ == "__main__":
    test_gui_screenshot()
