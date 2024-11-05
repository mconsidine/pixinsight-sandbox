import os
import json
import requests
import csv
import cv2
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget, QLabel, QPushButton, QHBoxLayout, QFileDialog, QGraphicsView,
    QGraphicsScene, QMessageBox, QInputDialog, QTreeWidget, QTreeWidgetItem, QCheckBox, QGraphicsTextItem, QDialog, QFormLayout, QSpinBox, QDialogButtonBox
)
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen, QColor, QTransform, QIcon
from PyQt5.QtCore import Qt, QRectF, QPointF
from PIL import Image
from astropy.io import fits
import tifffile as tiff
import sys
from scipy.interpolate import interp1d
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.coordinates import Angle
from astropy.wcs import WCS, Sip
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import time
from urllib.parse import quote
import urllib.parse
import webbrowser
from PyQt5.QtGui import QIcon
from math import sqrt
from astropy.utils.data import conf
conf.dataurl = 'file://C:/Users/Gaming/Desktop/Python Code/venv/Lib/site-packages/astroquery/simbad/data/'


# Access wrench_icon.png, adjusting for PyInstaller executable
if hasattr(sys, '_MEIPASS'):
    icon_path = os.path.join(sys._MEIPASS, 'wrench_icon.png')
    eye_icon_path = os.path.join(sys._MEIPASS, 'eye.png')
    disk_icon_path = os.path.join(sys._MEIPASS, 'disk.png')
    nuke_path = os.path.join(sys._MEIPASS, 'nuke.png')    
else:
    icon_path = 'wrench_icon.png'  # Path for running as a script
    eye_icon_path = 'eye.png'  # Path for running as a script
    disk_icon_path = 'disk.png'   
    nuke_path = 'nuke.png' 

# Constants for comoving radial distance calculation
H0 = 69.6  # Hubble constant in km/s/Mpc
WM = 0.286  # Omega(matter)
WV = 0.714  # Omega(vacuum)
c = 299792.458  # speed of light in km/s
Tyr = 977.8  # coefficient to convert 1/H into Gyr
Mpc_to_Gly = 3.262e-3  # Conversion from Mpc to Gly

otype_long_name_lookup = {
    "ev": "transient event",
    "Rad": "Radio-source",
    "mR": "metric Radio-source",
    "cm": "centimetric Radio-source",
    "mm": "millimetric Radio-source",
    "smm": "sub-millimetric source",
    "HI": "HI (21cm) source",
    "rB": "radio Burst",
    "Mas": "Maser",
    "IR": "Infra-Red source",
    "FIR": "Far-Infrared source",
    "MIR": "Mid-Infrared source",
    "NIR": "Near-Infrared source",
    "blu": "Blue object",
    "UV": "UV-emission source",
    "X": "X-ray source",
    "UX?": "Ultra-luminous X-ray candidate",
    "ULX": "Ultra-luminous X-ray source",
    "gam": "gamma-ray source",
    "gB": "gamma-ray Burst",
    "err": "Not an object (error, artefact, ...)",
    "grv": "Gravitational Source",
    "Lev": "(Micro)Lensing Event",
    "LS?": "Possible gravitational lens System",
    "Le?": "Possible gravitational lens",
    "LI?": "Possible gravitationally lensed image",
    "gLe": "Gravitational Lens",
    "gLS": "Gravitational Lens System (lens+images)",
    "GWE": "Gravitational Wave Event",
    "..?": "Candidate objects",
    "G?": "Possible Galaxy",
    "SC?": "Possible Supercluster of Galaxies",
    "C?G": "Possible Cluster of Galaxies",
    "Gr?": "Possible Group of Galaxies",
    "**?": "Physical Binary Candidate",
    "EB?": "Eclipsing Binary Candidate",
    "Sy?": "Symbiotic Star Candidate",
    "CV?": "Cataclysmic Binary Candidate",
    "No?": "Nova Candidate",
    "XB?": "X-ray binary Candidate",
    "LX?": "Low-Mass X-ray binary Candidate",
    "HX?": "High-Mass X-ray binary Candidate",
    "Pec?": "Possible Peculiar Star",
    "Y*?": "Young Stellar Object Candidate",
    "TT?": "T Tau star Candidate",
    "C*?": "Possible Carbon Star",
    "S*?": "Possible S Star",
    "OH?": "Possible Star with envelope of OH/IR type",
    "WR?": "Possible Wolf-Rayet Star",
    "Be?": "Possible Be Star",
    "Ae?": "Possible Herbig Ae/Be Star",
    "HB?": "Possible Horizontal Branch Star",
    "RR?": "Possible Star of RR Lyr type",
    "Ce?": "Possible Cepheid",
    "WV?": "Possible Variable Star of W Vir type",
    "RB?": "Possible Red Giant Branch star",
    "sg?": "Possible Supergiant star",
    "s?r": "Possible Red supergiant star",
    "s?y": "Possible Yellow supergiant star",
    "s?b": "Possible Blue supergiant star",
    "AB?": "Asymptotic Giant Branch Star candidate",
    "LP?": "Long Period Variable candidate",
    "Mi?": "Mira candidate",
    "pA?": "Post-AGB Star Candidate",
    "BS?": "Candidate blue Straggler Star",
    "HS?": "Hot subdwarf candidate",
    "WD?": "White Dwarf Candidate",
    "N*?": "Neutron Star Candidate",
    "BH?": "Black Hole Candidate",
    "SN?": "SuperNova Candidate",
    "LM?": "Low-mass star candidate",
    "BD?": "Brown Dwarf Candidate",
    "mul": "Composite object",
    "reg": "Region defined in the sky",
    "vid": "Underdense region of the Universe",
    "SCG": "Supercluster of Galaxies",
    "ClG": "Cluster of Galaxies",
    "GrG": "Group of Galaxies",
    "CGG": "Compact Group of Galaxies",
    "PaG": "Pair of Galaxies",
    "IG": "Interacting Galaxies",
    "C?*": "Possible (open) star cluster",
    "Gl?": "Possible Globular Cluster",
    "Cl*": "Cluster of Stars",
    "GlC": "Globular Cluster",
    "OpC": "Open (galactic) Cluster",
    "As*": "Association of Stars",
    "St*": "Stellar Stream",
    "MGr": "Moving Group",
    "**": "Double or multiple star",
    "EB*": "Eclipsing binary",
    "Al*": "Eclipsing binary of Algol type",
    "bL*": "Eclipsing binary of beta Lyr type",
    "WU*": "Eclipsing binary of W UMa type",
    "SB*": "Spectroscopic binary",
    "El*": "Ellipsoidal variable Star",
    "Sy*": "Symbiotic Star",
    "CV*": "Cataclysmic Variable Star",
    "DQ*": "CV DQ Her type (intermediate polar)",
    "AM*": "CV of AM Her type (polar)",
    "NL*": "Nova-like Star",
    "No*": "Nova",
    "DN*": "Dwarf Nova",
    "XB*": "X-ray Binary",
    "LXB": "Low Mass X-ray Binary",
    "HXB": "High Mass X-ray Binary",
    "ISM": "Interstellar matter",
    "PoC": "Part of Cloud",
    "PN?": "Possible Planetary Nebula",
    "CGb": "Cometary Globule",
    "bub": "Bubble",
    "EmO": "Emission Object",
    "Cld": "Cloud",
    "GNe": "Galactic Nebula",
    "DNe": "Dark Cloud (nebula)",
    "RNe": "Reflection Nebula",
    "MoC": "Molecular Cloud",
    "glb": "Globule (low-mass dark cloud)",
    "cor": "Dense core",
    "SFR": "Star forming region",
    "HVC": "High-velocity Cloud",
    "HII": "HII (ionized) region",
    "PN": "Planetary Nebula",
    "sh": "HI shell",
    "SR?": "SuperNova Remnant Candidate",
    "SNR": "SuperNova Remnant",
    "of?": "Outflow candidate",
    "out": "Outflow",
    "HH": "Herbig-Haro Object",
    "*": "Star",
    "V*?": "Star suspected of Variability",
    "Pe*": "Peculiar Star",
    "HB*": "Horizontal Branch Star",
    "Y*O": "Young Stellar Object",
    "Ae*": "Herbig Ae/Be star",
    "Em*": "Emission-line Star",
    "Be*": "Be Star",
    "BS*": "Blue Straggler Star",
    "RG*": "Red Giant Branch star",
    "AB*": "Asymptotic Giant Branch Star (He-burning)",
    "C*": "Carbon Star",
    "S*": "S Star",
    "sg*": "Evolved supergiant star",
    "s*r": "Red supergiant star",
    "s*y": "Yellow supergiant star",
    "s*b": "Blue supergiant star",
    "HS*": "Hot subdwarf",
    "pA*": "Post-AGB Star (proto-PN)",
    "WD*": "White Dwarf",
    "LM*": "Low-mass star (M<1solMass)",
    "BD*": "Brown Dwarf (M<0.08solMass)",
    "N*": "Confirmed Neutron Star",
    "OH*": "OH/IR star",
    "TT*": "T Tau-type Star",
    "WR*": "Wolf-Rayet Star",
    "PM*": "High proper-motion Star",
    "HV*": "High-velocity Star",
    "V*": "Variable Star",
    "Ir*": "Variable Star of irregular type",
    "Or*": "Variable Star of Orion Type",
    "Er*": "Eruptive variable Star",
    "RC*": "Variable Star of R CrB type",
    "RC?": "Variable Star of R CrB type candidate",
    "Ro*": "Rotationally variable Star",
    "a2*": "Variable Star of alpha2 CVn type",
    "Psr": "Pulsar",
    "BY*": "Variable of BY Dra type",
    "RS*": "Variable of RS CVn type",
    "Pu*": "Pulsating variable Star",
    "RR*": "Variable Star of RR Lyr type",
    "Ce*": "Cepheid variable Star",
    "dS*": "Variable Star of delta Sct type",
    "RV*": "Variable Star of RV Tau type",
    "WV*": "Variable Star of W Vir type",
    "bC*": "Variable Star of beta Cep type",
    "cC*": "Classical Cepheid (delta Cep type)",
    "gD*": "Variable Star of gamma Dor type",
    "SX*": "Variable Star of SX Phe type (subdwarf)",
    "LP*": "Long-period variable star",
    "Mi*": "Variable Star of Mira Cet type",
    "SN*": "SuperNova",
    "su*": "Sub-stellar object",
    "Pl?": "Extra-solar Planet Candidate",
    "Pl": "Extra-solar Confirmed Planet",
    "G": "Galaxy",
    "PoG": "Part of a Galaxy",
    "GiC": "Galaxy in Cluster of Galaxies",
    "BiC": "Brightest galaxy in a Cluster (BCG)",
    "GiG": "Galaxy in Group of Galaxies",
    "GiP": "Galaxy in Pair of Galaxies",
    "rG": "Radio Galaxy",
    "H2G": "HII Galaxy",
    "LSB": "Low Surface Brightness Galaxy",
    "AG?": "Possible Active Galaxy Nucleus",
    "Q?": "Possible Quasar",
    "Bz?": "Possible Blazar",
    "BL?": "Possible BL Lac",
    "EmG": "Emission-line galaxy",
    "SBG": "Starburst Galaxy",
    "bCG": "Blue compact Galaxy",
    "LeI": "Gravitationally Lensed Image",
    "LeG": "Gravitationally Lensed Image of a Galaxy",
    "LeQ": "Gravitationally Lensed Image of a Quasar",
    "AGN": "Active Galaxy Nucleus",
    "LIN": "LINER-type Active Galaxy Nucleus",
    "SyG": "Seyfert Galaxy",
    "Sy1": "Seyfert 1 Galaxy",
    "Sy2": "Seyfert 2 Galaxy",
    "Bla": "Blazar",
    "BLL": "BL Lac - type object",
    "OVV": "Optically Violently Variable object",
    "QSO": "Quasar"
}


# Configure Simbad to include the necessary fields, including redshift
Simbad.add_votable_fields('otype', 'otypes', 'diameter', 'z_value')
Simbad.ROW_LIMIT = 0  # Remove row limit for full results
Simbad.TIMEOUT = 60  # Increase timeout for long queries

# Astrometry.net API constants
ASTROMETRY_API_URL = "http://nova.astrometry.net/api/"
ASTROMETRY_API_KEY_FILE = "astrometry_api_key.txt"


def load_api_key():
    """Load the API key from a file, if it exists."""
    if os.path.exists(ASTROMETRY_API_KEY_FILE):
        with open(ASTROMETRY_API_KEY_FILE, 'r') as file:
            return file.read().strip()
    return None  # Return None if the file doesn't exist

def save_api_key(api_key):
    """Save the API key to a file."""
    with open(ASTROMETRY_API_KEY_FILE, 'w') as file:
        file.write(api_key)

def load_image(filename):
    try:
        file_extension = filename.lower().split('.')[-1]
        bit_depth = None
        is_mono = True
        original_header = None

        if file_extension in ['tif', 'tiff']:
            image = tiff.imread(filename)
            print(f"Loaded TIFF image with dtype: {image.dtype}")
            
            if image.dtype == np.uint16:
                image = image.astype(np.float32) / 65535.0
                bit_depth = "16-bit"
            elif image.dtype == np.uint32:
                image = image.astype(np.float32) / 4294967295.0
                bit_depth = "32-bit unsigned"
            else:
                image = image.astype(np.float32)
                bit_depth = "32-bit floating point"

            print(f"Final bit depth set to: {bit_depth}")

            # Check if the image has an alpha channel and remove it if necessary
            if image.shape[-1] == 4:
                print("Detected alpha channel in TIFF. Removing it.")
                image = image[:, :, :3]  # Keep only the first 3 channels (RGB)

            if len(image.shape) == 2:
                image = np.stack([image] * 3, axis=-1)
                is_mono = True

        elif file_extension in ['fits', 'fit']:
            with fits.open(filename) as hdul:
                image_data = hdul[0].data
                original_header = hdul[0].header  # Capture the FITS header

                # Determine the bit depth based on the data type in the FITS file
                if image_data.dtype == np.uint16:
                    bit_depth = "16-bit"
                    print("Identified 16-bit FITS image.")
                elif image_data.dtype == np.uint8:
                    bit_depth = "8-bit"
                    print("Identified 8-bit FITS image")
                elif image_data.dtype == np.float32:
                    bit_depth = "32-bit floating point"
                    print("Identified 32-bit floating point FITS image.")
                elif image_data.dtype == np.uint32:
                    bit_depth = "32-bit unsigned"
                    print("Identified 32-bit unsigned FITS image.")

                # Handle 3D FITS data (e.g., RGB or multi-layered data)
                if image_data.ndim == 3 and image_data.shape[0] == 3:
                    image = np.transpose(image_data, (1, 2, 0))  # Reorder to (height, width, channels)

                    if bit_depth == "16-bit":
                        image = image.astype(np.float32) / 65535.0  # Normalize to [0, 1] for 16-bit
                    elif bit_depth == "8-bit":
                        image = image.astype(np.float32) / 255.0    
                    elif bit_depth == "32-bit unsigned":
                        # Apply BSCALE and BZERO if present
                        bzero = original_header.get('BZERO', 0)
                        bscale = original_header.get('BSCALE', 1)
                        image = image.astype(np.float32) * bscale + bzero

                        # Normalize based on the actual data range
                        image_min, image_max = image.min(), image.max()
                        image = (image - image_min) / (image_max - image_min)
                        print(f"Image range after applying BZERO and BSCALE (3D case): min={image_min}, max={image_max}")

                    is_mono = False  # RGB data

                # Handle 2D FITS data (grayscale)
                elif image_data.ndim == 2:
                    if bit_depth == "16-bit":
                        image = image_data.astype(np.float32) / 65535.0  # Normalize to [0, 1] for 16-bit
                    elif bit_depth == "8-bit":
                        image = image_data.astype(np.float32) / 255.1    
                    elif bit_depth == "32-bit unsigned":
                        # Apply BSCALE and BZERO if present
                        bzero = original_header.get('BZERO', 0)
                        bscale = original_header.get('BSCALE', 1)
                        image = image_data.astype(np.float32) * bscale + bzero

                        # Normalize based on the actual data range
                        image_min, image_max = image.min(), image.max()
                        image = (image - image_min) / (image_max - image_min)
                        print(f"Image range after applying BZERO and BSCALE (2D case): min={image_min}, max={image_max}")

                    elif bit_depth == "32-bit floating point":
                        image = image_data  # No normalization needed for 32-bit float

                    is_mono = True
                    image = np.stack([image] * 3, axis=-1)  # Convert to 3-channel for consistency
                else:
                    raise ValueError("Unsupported FITS format!")

        else:
            # For PNG, JPEG, or other standard image formats
            image = np.array(Image.open(filename).convert('RGB')).astype(np.float32) / 255.0
            is_mono = False

        return image, original_header, bit_depth, is_mono

    except Exception as e:
        print(f"Error reading image {filename}: {e}")
        return None, None, None, None


class CustomGraphicsView(QGraphicsView):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.parent = parent
        self.setMouseTracking(True)  # Enable mouse tracking
        self.setDragMode(QGraphicsView.NoDrag)  # Disable default drag mode to avoid hand cursor
        self.setCursor(Qt.ArrowCursor)  # Set default cursor to arrow

        self.selected_object = None  # Initialize selected_object to None
        self.show_names = False 

        # Variables for drawing the circle
        self.circle_center = None
        self.circle_radius = 0
        self.drawing_circle = False  # Flag to check if we're currently drawing a circle
        self.dragging = False  # Flag to manage manual dragging


    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            if event.modifiers() == Qt.ShiftModifier:
                # Start drawing a circle for Shift+Click
                self.drawing_circle = True
                self.circle_center = self.mapToScene(event.pos())
                self.circle_radius = 0
                self.parent.status_label.setText("Drawing circle: Shift + Drag")
                self.update_circle()
            else:
                # Detect if an object circle was clicked without Shift
                scene_pos = self.mapToScene(event.pos())
                clicked_object = self.get_object_at_position(scene_pos)
                
                if clicked_object:
                    # Select the clicked object and redraw
                    self.select_object(clicked_object)
                    self.draw_query_results()
                    
                    # Highlight the corresponding row in the TreeWidget
                    for i in range(self.parent.results_tree.topLevelItemCount()):
                        item = self.parent.results_tree.topLevelItem(i)
                        if item.text(2) == clicked_object["name"]:   # Assuming third element is 'Name'
                            self.parent.results_tree.setCurrentItem(item)
                            break
                else:
                    # Start manual dragging if no Shift is held
                    self.dragging = True
                    self.setCursor(Qt.ClosedHandCursor)  # Use closed hand cursor to indicate dragging
                    self.drag_start_pos = event.pos()  # Store starting position
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event):
        """Handle double-click event on an object in the main image to open SIMBAD or NED URL based on source."""
        scene_pos = self.mapToScene(event.pos())
        clicked_object = self.get_object_at_position(scene_pos)

        if clicked_object:
            object_name = clicked_object.get("name")  # Access 'name' key from the dictionary
            ra = float(clicked_object.get("ra"))  # Ensure RA is a float for precision
            dec = float(clicked_object.get("dec"))  # Ensure Dec is a float for precision
            source = clicked_object.get("source", "Simbad")  # Default to "Simbad" if source not specified

            if source == "Simbad" and object_name:
                # Open Simbad URL with encoded object name
                encoded_name = quote(object_name)
                url = f"https://simbad.cds.unistra.fr/simbad/sim-basic?Ident={encoded_name}&submit=SIMBAD+search"
                webbrowser.open(url)
            elif source == "Vizier":
                # Format the NED search URL with proper RA, Dec, and radius
                radius = 5 / 60  # Radius in arcminutes (5 arcseconds)
                dec_sign = "%2B" if dec >= 0 else "-"  # Determine sign for declination
                ned_url = (
                    f"http://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search"
                    f"&ra={ra:.6f}d&dec={dec_sign}{abs(dec):.6f}d&radius={radius:.3f}"
                    "&in_csys=Equatorial&in_equinox=J2000.0"
                )
                webbrowser.open(ned_url)
        else:
            super().mouseDoubleClickEvent(event)

    def mouseMoveEvent(self, event):
        if self.drawing_circle:
            # Update the circle radius as the mouse moves
            scene_pos = self.mapToScene(event.pos())
            self.circle_radius = np.sqrt(
                (scene_pos.x() - self.circle_center.x()) ** 2 +
                (scene_pos.y() - self.circle_center.y()) ** 2
            )
            self.update_circle()
        elif self.dragging:
            # Handle manual dragging by scrolling the view
            delta = event.pos() - self.drag_start_pos
            self.horizontalScrollBar().setValue(self.horizontalScrollBar().value() - delta.x())
            self.verticalScrollBar().setValue(self.verticalScrollBar().value() - delta.y())
            self.drag_start_pos = event.pos()  # Update starting position for smooth dragging
        else:
            # Update RA/Dec display as the cursor moves
            self.parent.update_ra_dec_from_mouse(event)
        super().mouseMoveEvent(event)
            

    def mouseReleaseEvent(self, event):
        if self.drawing_circle and event.button() == Qt.LeftButton:
            # Stop drawing the circle
            self.drawing_circle = False
            self.parent.circle_center = self.circle_center
            self.parent.circle_radius = self.circle_radius

            # Calculate RA/Dec for the circle center
            ra, dec = self.parent.calculate_ra_dec_from_pixel(self.circle_center.x(), self.circle_center.y())
            if ra is not None and dec is not None:
                self.parent.ra_label.setText(f"RA: {self.parent.convert_ra_to_hms(ra)}")
                self.parent.dec_label.setText(f"Dec: {self.parent.convert_dec_to_dms(dec)}")

                if self.parent.pixscale:
                    radius_arcmin = self.circle_radius * self.parent.pixscale / 60.0
                    self.parent.status_label.setText(
                        f"Circle set at center RA={ra:.6f}, Dec={dec:.6f}, radius={radius_arcmin:.2f} arcmin"
                    )
                else:
                    self.parent.status_label.setText("Pixscale not available for radius calculation.")
            else:
                self.parent.status_label.setText("Unable to determine RA/Dec due to missing WCS.")

            # Update circle data and redraw
            self.parent.update_circle_data()
            self.update_circle()

        # Stop manual dragging and reset cursor to arrow
        self.dragging = False
        self.setCursor(Qt.ArrowCursor)
        
        # Update the mini preview to reflect any changes
        self.update_mini_preview()

        super().mouseReleaseEvent(event)

    def wheelEvent(self, event):
        """Handle zoom in and out with the mouse wheel."""
        if event.angleDelta().y() > 0:
            self.parent.zoom_in()
        else:
            self.parent.zoom_out()        

    def update_circle(self):
        if self.parent.main_image:
            self.parent.main_scene.clear()
            self.parent.main_scene.addPixmap(self.parent.main_image)
            
            pen_circle = QPen(QColor(255, 0, 0), 2)
            self.parent.main_scene.addEllipse(
                int(self.circle_center.x() - self.circle_radius),
                int(self.circle_center.y() - self.circle_radius),
                int(self.circle_radius * 2),
                int(self.circle_radius * 2),
                pen_circle
            )
            self.update_mini_preview()

    def scrollContentsBy(self, dx, dy):
        """Called whenever the main preview scrolls, ensuring the green box updates in the mini preview."""
        super().scrollContentsBy(dx, dy)
        self.parent.update_green_box()

    def update_mini_preview(self):
        """Update the mini preview with the current view rectangle and any additional mirrored elements."""
        if self.parent.main_image:
            # Scale the main image to fit in the mini preview
            mini_pixmap = self.parent.main_image.scaled(
                self.parent.mini_preview.size(),
                Qt.KeepAspectRatio,
                Qt.SmoothTransformation
            )
            mini_painter = QPainter(mini_pixmap)

            try:
                # Define scale factors based on main image dimensions
                if self.parent.main_image.width() > 0 and self.parent.main_image.height() > 0:
                    scale_factor_x = mini_pixmap.width() / self.parent.main_image.width()
                    scale_factor_y = mini_pixmap.height() / self.parent.main_image.height()

                    # Draw the search circle if it's defined
                    if self.circle_center is not None and self.circle_radius > 0:
                        pen_circle = QPen(QColor(255, 0, 0), 2)
                        mini_painter.setPen(pen_circle)
                        mini_painter.drawEllipse(
                            int(self.circle_center.x() * scale_factor_x - self.circle_radius * scale_factor_x),
                            int(self.circle_center.y() * scale_factor_y - self.circle_radius * scale_factor_y),
                            int(self.circle_radius * 2 * scale_factor_x),
                            int(self.circle_radius * 2 * scale_factor_y)
                        )

                    # Draw the green box representing the current view
                    mini_painter.setPen(QPen(QColor(0, 255, 0), 2))
                    view_rect = self.parent.main_preview.mapToScene(
                        self.parent.main_preview.viewport().rect()
                    ).boundingRect()
                    mini_painter.drawRect(
                        int(view_rect.x() * scale_factor_x),
                        int(view_rect.y() * scale_factor_y),
                        int(view_rect.width() * scale_factor_x),
                        int(view_rect.height() * scale_factor_y)
                    )

                    # Draw red dots for each result instead of green circles
                    mini_painter.setPen(QPen(QColor(255, 0, 0), 4))  # Red color with width 4 for dots
                    for obj in self.parent.results:
                        ra, dec = obj['ra'], obj['dec']
                        x, y = self.parent.calculate_pixel_from_ra_dec(ra, dec)
                        if x is not None and y is not None:
                            mini_painter.drawPoint(
                                int(x * scale_factor_x),
                                int(y * scale_factor_y)
                            )
            finally:
                mini_painter.end()  # Ensure QPainter is properly ended

            self.parent.mini_preview.setPixmap(mini_pixmap)

    def draw_query_results(self):
        """Draw query results with or without names based on the show_names setting."""
        if self.parent.main_image:
            # Clear the main scene and re-add the main image
            self.parent.main_scene.clear()
            self.parent.main_scene.addPixmap(self.parent.main_image)
            
            # Ensure the search circle is drawn
            self.update_circle()

            for obj in self.parent.results:
                ra, dec, name = obj["ra"], obj["dec"], obj["name"]
                x, y = self.parent.calculate_pixel_from_ra_dec(ra, dec)
                if x is not None and y is not None:
                    # Use green color if this is the selected object, otherwise use red
                    pen_color = QColor(0, 255, 0) if obj == self.selected_object else QColor(255, 0, 0)
                    pen = QPen(pen_color, 2)
                    self.parent.main_scene.addEllipse(int(x - 5), int(y - 5), 10, 10, pen)

                    # Conditionally draw names if the checkbox is checked
                    if self.parent.show_names:
                        text_item = QGraphicsTextItem(name)
                        text_item.setPos(x + 10, y + 10)
                        text_item.setDefaultTextColor(QColor(255, 255, 0))  # Set text color to white
                        self.parent.main_scene.addItem(text_item)          

    def clear_query_results(self):
        """Clear markers from the main image."""
        self.parent.main_scene.clear()
        if self.parent.main_image:
            self.parent.main_scene.addPixmap(self.parent.main_image)
            self.parent.update_circle_data()                        

    def set_query_results(self, results):
        """Set the query results and redraw."""
        self.parent.results = results  # Store results as dictionaries
        self.draw_query_results()

    def get_object_at_position(self, pos):
        """Find the object at the given position in the main preview."""
        for obj in self.parent.results:
            ra, dec = obj["ra"], obj["dec"]
            x, y = self.parent.calculate_pixel_from_ra_dec(ra, dec)
            if x is not None and y is not None:
                if abs(pos.x() - x) <= 5 and abs(pos.y() - y) <= 5:
                    return obj
        return None


    def select_object(self, selected_obj):
        """Select or deselect the specified object and update visuals."""
        self.selected_object = selected_obj if self.selected_object != selected_obj else None
        self.draw_query_results()  # Redraw to reflect selection

        # Update the TreeWidget selection in MainWindow
        for i in range(self.parent.results_tree.topLevelItemCount()):
            item = self.parent.results_tree.topLevelItem(i)
            if item.text(2) == selected_obj["name"]:  # Assuming 'name' is the unique identifier
                self.parent.results_tree.setCurrentItem(item if self.selected_object else None)
                break

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("What's In My Image - Standalone")
        self.setGeometry(100, 100, 1200, 800)
        self.circle_center = None
        self.circle_radius = 0    
        self.show_names = False  # Boolean to toggle showing names on the main image
        self.max_results = 100  # Default maximum number of query results                 

        main_layout = QHBoxLayout()

        # Left Column Layout
        left_panel = QVBoxLayout()
        
        # Load button
        self.load_button = QPushButton("Load Image")
        self.load_button.clicked.connect(self.open_image)
        left_panel.addWidget(self.load_button)

        # RA/Dec Labels
        ra_dec_layout = QHBoxLayout()
        self.ra_label = QLabel("RA: N/A")
        self.dec_label = QLabel("Dec: N/A")
        ra_dec_layout.addWidget(self.ra_label)
        ra_dec_layout.addWidget(self.dec_label)
        left_panel.addLayout(ra_dec_layout)

        # Query Simbad button
        self.query_button = QPushButton("Query Simbad")
        left_panel.addWidget(self.query_button)
        self.query_button.clicked.connect(lambda: self.query_simbad(self.get_defined_radius()))


        # Create a horizontal layout for the show names checkbox and clear results button
        show_clear_layout = QHBoxLayout()

        # Create the Show Object Names checkbox
        self.show_names_checkbox = QCheckBox("Show Object Names")
        self.show_names_checkbox.stateChanged.connect(self.toggle_object_names)  # Connect to a function to toggle names
        show_clear_layout.addWidget(self.show_names_checkbox)

        # Create the Clear Results button
        self.clear_results_button = QPushButton("Clear Results")
        self.clear_results_button.clicked.connect(self.clear_search_results)  # Connect to a function to clear results
        show_clear_layout.addWidget(self.clear_results_button)

        # Add this horizontal layout to the left panel layout (or wherever you want it to appear)
        left_panel.addLayout(show_clear_layout)   

        # Create a horizontal layout for the two buttons
        button_layout = QHBoxLayout()

        # Show Visible Objects Only button
        self.toggle_visible_objects_button = QPushButton("Show Visible Objects Only")
        self.toggle_visible_objects_button.setCheckable(True)  # Toggle button state
        self.toggle_visible_objects_button.setIcon(QIcon(eye_icon_path))
        self.toggle_visible_objects_button.clicked.connect(self.filter_visible_objects)
        self.toggle_visible_objects_button.setToolTip("Toggle the visibility of objects based on brightness.")
        button_layout.addWidget(self.toggle_visible_objects_button)

        # Save CSV button
        self.save_csv_button = QPushButton("Save CSV")
        self.save_csv_button.setIcon(QIcon(disk_icon_path))
        self.save_csv_button.clicked.connect(self.save_results_as_csv)
        button_layout.addWidget(self.save_csv_button)

        # Add the button layout to the left panel or main layout
        left_panel.addLayout(button_layout)  

        # Advanced Search Button
        self.advanced_search_button = QPushButton("Advanced Search")
        self.advanced_search_button.setCheckable(True)
        self.advanced_search_button.clicked.connect(self.toggle_advanced_search)
        left_panel.addWidget(self.advanced_search_button)

        # Advanced Search Panel (initially hidden)
        self.advanced_search_panel = QVBoxLayout()
        self.advanced_search_panel_widget = QWidget()
        self.advanced_search_panel_widget.setLayout(self.advanced_search_panel)
        self.advanced_search_panel_widget.setFixedWidth(300)
        self.advanced_search_panel_widget.setVisible(False)  # Hide initially        

        # Status label
        self.status_label = QLabel("Status: Ready")
        left_panel.addWidget(self.status_label)

        # Settings (wrench icon) button
        self.settings_button = QPushButton()
        self.settings_button.setIcon(QIcon(icon_path))  # Ensure you have this icon in the project directory
        self.settings_button.clicked.connect(self.open_settings_dialog)
        left_panel.addWidget(self.settings_button)        

        # Mini Preview
        self.mini_preview = QLabel("Mini Preview")
        self.mini_preview.setFixedSize(300, 300)
        self.mini_preview.mousePressEvent = self.on_mini_preview_press
        self.mini_preview.mouseMoveEvent = self.on_mini_preview_drag
        self.mini_preview.mouseReleaseEvent = self.on_mini_preview_release
        left_panel.addWidget(self.mini_preview)

        # Right Column Layout
        right_panel = QVBoxLayout()

        # Zoom buttons above the main preview
        zoom_controls_layout = QHBoxLayout()
        self.zoom_in_button = QPushButton("Zoom In")
        self.zoom_in_button.clicked.connect(self.zoom_in)
        self.zoom_out_button = QPushButton("Zoom Out")
        self.zoom_out_button.clicked.connect(self.zoom_out)
        zoom_controls_layout.addWidget(self.zoom_in_button)
        zoom_controls_layout.addWidget(self.zoom_out_button)
        right_panel.addLayout(zoom_controls_layout)        

        # Main Preview
        self.main_preview = CustomGraphicsView(self)
        self.main_scene = QGraphicsScene(self.main_preview)
        self.main_preview.setScene(self.main_scene)
        self.main_preview.setRenderHint(QPainter.Antialiasing)
        self.main_preview.setTransformationAnchor(QGraphicsView.AnchorUnderMouse)
        right_panel.addWidget(self.main_preview)

        # Connect scroll events to update the green box in the mini preview
        self.main_preview.verticalScrollBar().valueChanged.connect(self.main_preview.update_mini_preview)
        self.main_preview.horizontalScrollBar().valueChanged.connect(self.main_preview.update_mini_preview)

        # Tree Widget for results
        self.results_tree = QTreeWidget()
        self.results_tree.setHeaderLabels(["RA", "Dec", "Name", "Diameter", "Type", "Long Type", "Redshift", "Comoving Radial Distance (GLy)"])
        self.results_tree.setFixedHeight(150)
        self.results_tree.itemClicked.connect(self.on_tree_item_clicked)
        self.results_tree.itemDoubleClicked.connect(self.on_tree_item_double_clicked) 
        self.results_tree.setSortingEnabled(True)  # <-- Add this line
        right_panel.addWidget(self.results_tree)

        # Advanced Search Panel
        self.advanced_param_label = QLabel("Advanced Search Parameters")
        self.advanced_search_panel.addWidget(self.advanced_param_label)

        # TreeWidget for object types
        self.object_tree = QTreeWidget()
        self.object_tree.setHeaderLabels(["Object Type", "Description"])
        self.object_tree.setColumnWidth(0, 150)
        self.object_tree.setSortingEnabled(True)

        # Populate the TreeWidget with object types from otype_long_name_lookup
        for obj_type, description in otype_long_name_lookup.items():
            item = QTreeWidgetItem([obj_type, description])
            item.setCheckState(0, Qt.Checked)  # Start with all items unchecked
            self.object_tree.addTopLevelItem(item)

        self.advanced_search_panel.addWidget(self.object_tree)

        # Buttons for toggling selections
        toggle_buttons_layout = QHBoxLayout()

        # Toggle All Button
        self.toggle_all_button = QPushButton("Toggle All")
        self.toggle_all_button.clicked.connect(self.toggle_all_items)
        toggle_buttons_layout.addWidget(self.toggle_all_button)

        # Toggle Stars Button
        self.toggle_stars_button = QPushButton("Toggle Stars")
        self.toggle_stars_button.clicked.connect(self.toggle_star_items)
        toggle_buttons_layout.addWidget(self.toggle_stars_button)

        # Toggle Galaxies Button
        self.toggle_galaxies_button = QPushButton("Toggle Galaxies")
        self.toggle_galaxies_button.clicked.connect(self.toggle_galaxy_items)
        toggle_buttons_layout.addWidget(self.toggle_galaxies_button)

        # Add toggle buttons to the advanced search layout
        self.advanced_search_panel.addLayout(toggle_buttons_layout)    

        # Add Simbad Search buttons below the toggle buttons
        search_button_layout = QHBoxLayout()

        self.simbad_defined_region_button = QPushButton("Search Defined Region")
        self.simbad_defined_region_button.clicked.connect(self.search_defined_region)
        search_button_layout.addWidget(self.simbad_defined_region_button)

        self.simbad_entire_image_button = QPushButton("Search Entire Image")
        self.simbad_entire_image_button.clicked.connect(self.search_entire_image)
        search_button_layout.addWidget(self.simbad_entire_image_button)

        self.advanced_search_panel.addLayout(search_button_layout)

        # Adding the "Deep Vizier Search" button below the other search buttons
        self.deep_vizier_button = QPushButton("Caution - Deep Vizier Search")
        self.deep_vizier_button.setIcon(QIcon(nuke_path))  # Assuming `nuke_path` is the correct path for the icon
        self.deep_vizier_button.setToolTip("Perform a deep search with Vizier. Caution: May return large datasets.")

        # Connect the button to a placeholder method for the deep Vizier search
        self.deep_vizier_button.clicked.connect(self.perform_deep_vizier_search)

        # Add the Deep Vizier button to the advanced search layout
        self.advanced_search_panel.addWidget(self.deep_vizier_button)
                        

        # Combine left and right panels
        main_layout.addLayout(left_panel)
        main_layout.addLayout(right_panel)
        main_layout.addWidget(self.advanced_search_panel_widget)
        
        container = QWidget()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        self.image_path = None
        self.zoom_level = 1.0
        self.main_image = None
        self.green_box = None
        self.dragging = False
        self.center_ra = None
        self.center_dec = None
        self.pixscale = None
        self.orientation = None
        self.parity = None  
        self.circle_center = None
        self.circle_radius = 0  
        self.results = []
        self.wcs = None  # Initialize WCS to None

    def get_selected_object_types(self):
        """Return a list of selected object types from the tree widget."""
        selected_types = []
        for i in range(self.object_tree.topLevelItemCount()):
            item = self.object_tree.topLevelItem(i)
            if item.checkState(0) == Qt.Checked:
                selected_types.append(item.text(0))  # Add the object type
        return selected_types
    
    def search_defined_region(self):
        """Perform a Simbad search for the defined region and filter by selected object types."""
        selected_types = self.get_selected_object_types()
        if not selected_types:
            QMessageBox.warning(self, "No Object Types Selected", "Please select at least one object type.")
            return

        # Calculate the radius in degrees for the defined region (circle radius)
        radius_deg = self.get_defined_radius()

        # Perform the Simbad search in the defined region with the calculated radius
        self.query_simbad(radius_deg)


    def search_entire_image(self):
        """Search the entire image using Simbad with selected object types."""
        selected_types = self.get_selected_object_types()  # Get selected types from the advanced search panel
        if not selected_types:
            QMessageBox.warning(self, "No Object Types Selected", "Please select at least one object type.")
            return

        # Calculate radius as the distance from the image center to a corner
        width, height = self.main_image.width(), self.main_image.height()
        center_x, center_y = width / 2, height / 2
        corner_x, corner_y = width, height  # Bottom-right corner
        # Calculate distance in pixels from center to corner
        radius_px = np.sqrt((corner_x - center_x) ** 2 + (corner_y - center_y) ** 2)
        # Convert radius from pixels to degrees
        radius_deg = float((radius_px * self.pixscale) / 3600.0)

        # Perform the query with the calculated radius
        self.query_simbad(radius_deg, max_results=25000)




    def toggle_advanced_search(self):
        """Toggle visibility of the advanced search panel."""
        self.advanced_search_panel.setVisible(not self.advanced_search_panel.isVisible())

    def toggle_all_items(self):
        """Toggle selection for all items in the object tree."""
        for i in range(self.object_tree.topLevelItemCount()):
            item = self.object_tree.topLevelItem(i)
            new_state = Qt.Checked if item.checkState(0) == Qt.Unchecked else Qt.Unchecked
            item.setCheckState(0, new_state)

    def toggle_star_items(self):
        """Toggle selection for items related to stars."""
        star_keywords = ["star", "Eclipsing binary of W UMa type", "Spectroscopic binary",
                         "Variable of RS CVn type", "Mira candidate", "Long Period Variable candidate",
                         "Hot subdwarf", "Eclipsing Binary Candidate", "Eclipsing binary", 
                         "Cataclysmic Binary Candidate", "Possible Cepheid", "White Dwarf", 
                         "White Dwarf Candidate"]
        for i in range(self.object_tree.topLevelItemCount()):
            item = self.object_tree.topLevelItem(i)
            description = item.text(1).lower()
            object_type = item.text(0)
            if any(keyword.lower() in description for keyword in star_keywords) or "*" in object_type:
                new_state = Qt.Checked if item.checkState(0) == Qt.Unchecked else Qt.Unchecked
                item.setCheckState(0, new_state)

    def toggle_galaxy_items(self):
        """Toggle selection for items related to galaxies."""
        for i in range(self.object_tree.topLevelItemCount()):
            item = self.object_tree.topLevelItem(i)
            description = item.text(1).lower()
            if "galaxy" in description or "galaxies" in description:
                new_state = Qt.Checked if item.checkState(0) == Qt.Unchecked else Qt.Unchecked
                item.setCheckState(0, new_state)


    def toggle_advanced_search(self):
        """Toggle the visibility of the advanced search panel."""
        is_visible = self.advanced_search_panel_widget.isVisible()
        self.advanced_search_panel_widget.setVisible(not is_visible)

    def save_results_as_csv(self):
        """Save the results from the TreeWidget as a CSV file."""
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if path:
            with open(path, mode='w', newline='') as file:
                writer = csv.writer(file)
                # Write header
                writer.writerow(["RA", "Dec", "Name", "Diameter", "Type", "Long Type", "Redshift", "Comoving Radial Distance (GLy)"])

                # Write data from TreeWidget
                for i in range(self.results_tree.topLevelItemCount()):
                    item = self.results_tree.topLevelItem(i)
                    row_data = [item.text(column) for column in range(self.results_tree.columnCount())]
                    writer.writerow(row_data)

            QMessageBox.information(self, "CSV Saved", f"Results successfully saved to {path}")        

    def filter_visible_objects(self):
        """Filter objects based on visibility threshold."""
        if not self.main_image:  # Ensure there's an image loaded
            QMessageBox.warning(self, "No Image", "Please load an image first.")
            return

        n = 0.2  # Threshold multiplier, adjust as needed
        median, std_dev = self.calculate_image_statistics(self.main_image)

        # Remove objects below threshold from results
        filtered_results = []
        for obj in self.results:
            if self.is_marker_visible(obj, median, std_dev, n):
                filtered_results.append(obj)

        # Update the results and redraw the markers
        self.results = filtered_results
        self.update_results_tree()
        self.main_preview.draw_query_results()

    def calculate_image_statistics(self, image):
        """Calculate median and standard deviation for a grayscale image efficiently using OpenCV."""
        
        # Convert QPixmap to QImage if necessary
        qimage = image.toImage()

        # Convert QImage to a format compatible with OpenCV
        width = qimage.width()
        height = qimage.height()
        ptr = qimage.bits()
        ptr.setsize(height * width * 4)  # 4 channels (RGBA)
        img_array = np.array(ptr).reshape(height, width, 4)  # Convert to RGBA array

        # Convert to grayscale for analysis
        gray_image = cv2.cvtColor(img_array, cv2.COLOR_RGBA2GRAY)

        # Calculate median and standard deviation
        median = np.median(gray_image)
        _, std_dev = cv2.meanStdDev(gray_image)

        return median, std_dev[0][0]  # std_dev returns a 2D array, so we extract the single value
    
    def is_marker_visible(self, marker, median, std_dev, n):
        """Check if the marker's brightness is above the threshold."""
        threshold = median + n * std_dev
        check_size = 8  # Define a 4x4 region around the marker

        # Convert QPixmap to QImage to access pixel colors
        image = self.main_image.toImage()

        # Get marker coordinates in pixel space
        ra, dec = marker.get('ra'), marker.get('dec')
        if ra is not None and dec is not None:
            x, y = self.calculate_pixel_from_ra_dec(ra, dec)
            if x is None or y is None:
                return False  # Skip marker if it can't be converted to pixels
        else:
            return False

        # Calculate brightness in a 4x4 region around marker coordinates
        brightness_values = []
        for dx in range(-check_size // 2, check_size // 2):
            for dy in range(-check_size // 2, check_size // 2):
                px = x + dx
                py = y + dy
                if 0 <= px < image.width() and 0 <= py < image.height():
                    color = image.pixelColor(px, py)  # Get color from QImage
                    brightness = color.value() if color.isValid() else 0  # Adjust for grayscale
                    brightness_values.append(brightness)

        if brightness_values:
            average_brightness = sum(brightness_values) / len(brightness_values)
            return average_brightness > threshold
        else:
            return False



    def update_results_tree(self):
        """Refresh the TreeWidget to reflect current results."""
        self.results_tree.clear()
        for obj in self.results:
            item = QTreeWidgetItem([
                str(obj['ra']),
                str(obj['dec']),
                obj['name'],
                str(obj['diameter']),
                obj['short_type'],
                obj['long_type'],
                str(obj['redshift']),
                str(obj['comoving_distance'])
            ])
            self.results_tree.addTopLevelItem(item)

    def toggle_object_names(self, state):
        """Toggle the visibility of object names based on the checkbox state."""
        self.show_names = state == Qt.Checked
        self.main_preview.draw_query_results()  # Redraw to apply the change


    # Function to clear search results and remove markers
    def clear_search_results(self):
        """Clear the search results and remove all markers."""
        self.results_tree.clear()        # Clear the results from the tree
        self.results = []                # Clear the results list
        self.main_preview.results = []   # Clear results from the main preview
        self.main_preview.selected_object = None
        self.main_preview.draw_query_results()  # Redraw the main image without markers
        self.status_label.setText("Results cleared.")

    def on_tree_item_clicked(self, item):
        """Handle item click in the TreeWidget to highlight the associated object."""
        object_name = item.text(2)

        # Find the object in results
        selected_object = next(
            (obj for obj in self.results if obj.get("name") == object_name), None
        )

        if selected_object:
            # Set the selected object in the CustomGraphicsView and redraw
            self.main_preview.select_object(selected_object)
            self.main_preview.draw_query_results()

    def on_tree_item_double_clicked(self, item):
        """Handle double-click event on a TreeWidget item to open SIMBAD or NED URL based on source."""
        object_name = item.text(2)  # Assuming 'Name' is in the third column
        ra = float(item.text(0).strip())  # Assuming RA is in the first column
        dec = float(item.text(1).strip())  # Assuming Dec is in the second column
        
        # Retrieve the entry directly from self.query_results
        entry = next((result for result in self.query_results if float(result['ra']) == ra and float(result['dec']) == dec), None)
        source = entry.get('source', 'Simbad') if entry else 'Simbad'  # Default to "Simbad" if entry not found

        if source == "Simbad" and object_name:
            # Open Simbad URL with encoded object name
            encoded_name = quote(object_name)
            simbad_url = f"https://simbad.cds.unistra.fr/simbad/sim-basic?Ident={encoded_name}&submit=SIMBAD+search"
            webbrowser.open(simbad_url)
        elif source == "Vizier":
            # Format the NED search URL with proper RA, Dec, and radius
            radius = 5 / 60  # Radius in arcminutes (5 arcseconds)
            dec_sign = "%2B" if dec >= 0 else "-"  # Determine sign for declination
            ned_url = f"http://ned.ipac.caltech.edu/conesearch?search_type=Near%20Position%20Search&ra={ra:.6f}d&dec={dec_sign}{abs(dec):.6f}d&radius={radius:.3f}&in_csys=Equatorial&in_equinox=J2000.0"
            webbrowser.open(ned_url)

    

    def open_image(self):
        self.image_path, _ = QFileDialog.getOpenFileName(self, "Open Image", "", "Images (*.png *.jpg *.jpeg *.tif *.fits)")
        if self.image_path:
            img_array, original_header, bit_depth, is_mono = load_image(self.image_path)
            if img_array is not None:
                # Prepare image for display
                img = (img_array * 255).astype(np.uint8)
                height, width, _ = img.shape
                bytes_per_line = 3 * width
                qimg = QImage(img.tobytes(), width, height, bytes_per_line, QImage.Format_RGB888)
                pixmap = QPixmap.fromImage(qimg)

                self.main_image = pixmap
                scaled_pixmap = pixmap.scaled(self.mini_preview.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
                self.mini_preview.setPixmap(scaled_pixmap)

                self.main_scene.clear()
                self.main_scene.addPixmap(pixmap)
                self.main_preview.setSceneRect(QRectF(pixmap.rect()))
                self.zoom_level = 1.0
                self.main_preview.resetTransform()
                self.main_preview.centerOn(self.main_scene.sceneRect().center())
                self.update_green_box()

                # Initialize WCS from FITS header if it is a FITS file
                if self.image_path.lower().endswith(('.fits', '.fit')):
                    with fits.open(self.image_path) as hdul:
                        self.header = hdul[0].header
                        
                        try:
                            # Use only the first two dimensions for WCS
                            self.wcs = WCS(self.header, naxis=2)
                            
                            # Calculate and set pixel scale
                            pixel_scale_matrix = self.wcs.pixel_scale_matrix
                            self.pixscale = np.sqrt(pixel_scale_matrix[0, 0]**2 + pixel_scale_matrix[1, 0]**2) * 3600  # arcsec/pixel
                            self.center_ra, self.center_dec = self.wcs.wcs.crval
                            self.print_corner_coordinates()
                            
                            # Display WCS information
                            print(f"WCS data loaded from FITS header: RA={self.center_ra}, Dec={self.center_dec}, "
                                f"Pixel Scale={self.pixscale} arcsec/px")
                            
                            
                        except ValueError as e:
                            print("Error initializing WCS:", e)
                            QMessageBox.warning(self, "WCS Error", "Failed to load WCS data from FITS header.")
                else:
                    # For non-FITS images (e.g., JPEG, PNG), prompt directly for a blind solve
                    self.prompt_blind_solve()


    def print_corner_coordinates(self):
        """Print the RA/Dec coordinates of the four corners of the image for debugging purposes."""
        if not hasattr(self, 'wcs'):
            print("WCS data is incomplete, cannot calculate corner coordinates.")
            return

        width = self.main_image.width()
        height = self.main_image.height()

        # Define the corner coordinates
        corners = {
            "Top-Left": (0, 0),
            "Top-Right": (width, 0),
            "Bottom-Left": (0, height),
            "Bottom-Right": (width, height)
        }

        print("Corner RA/Dec coordinates:")
        for corner_name, (x, y) in corners.items():
            ra, dec = self.calculate_ra_dec_from_pixel(x, y)
            ra_hms = self.convert_ra_to_hms(ra)
            dec_dms = self.convert_dec_to_dms(dec)
            print(f"{corner_name}: RA={ra_hms}, Dec={dec_dms}")

    def calculate_ra_dec_from_pixel(self, x, y):
        """Convert pixel coordinates (x, y) to RA/Dec using Astropy WCS."""
        if not hasattr(self, 'wcs'):
            print("WCS not initialized.")
            return None, None

        # Convert pixel coordinates to sky coordinates
        ra, dec = self.wcs.all_pix2world(x, y, 0)

        return ra, dec
                        


    def update_ra_dec_from_mouse(self, event):
        """Update RA and Dec based on mouse position over the main preview."""
        if self.main_image and self.wcs:
            pos = self.main_preview.mapToScene(event.pos())
            x, y = int(pos.x()), int(pos.y())
            
            if 0 <= x < self.main_image.width() and 0 <= y < self.main_image.height():
                # Convert using the saved RA/Dec if it’s already calculated
                ra, dec = self.calculate_ra_dec_from_pixel(x, y)
                ra_hms = self.convert_ra_to_hms(ra)
                dec_dms = self.convert_dec_to_dms(dec)

                # Update labels
                self.ra_label.setText(f"RA: {ra_hms}")
                self.dec_label.setText(f"Dec: {dec_dms}")
        else:
            # If WCS data is missing, update labels to indicate RA/Dec is unavailable
            self.ra_label.setText("RA: N/A")
            self.dec_label.setText("Dec: N/A")


    def convert_ra_to_hms(self, ra_deg):
        """Convert Right Ascension in degrees to Hours:Minutes:Seconds format."""
        ra_hours = ra_deg / 15.0  # Convert degrees to hours
        hours = int(ra_hours)
        minutes = int((ra_hours - hours) * 60)
        seconds = (ra_hours - hours - minutes / 60.0) * 3600
        return f"{hours:02d}h{minutes:02d}m{seconds:05.2f}s"

    def convert_dec_to_dms(self, dec_deg):
        """Convert Declination in degrees to Degrees:Minutes:Seconds format."""
        sign = "-" if dec_deg < 0 else "+"
        dec_deg = abs(dec_deg)
        degrees = int(dec_deg)
        minutes = int((dec_deg - degrees) * 60)
        seconds = (dec_deg - degrees - minutes / 60.0) * 3600
        degree_symbol = "\u00B0"
        return f"{sign}{degrees:02d}{degree_symbol}{minutes:02d}m{seconds:05.2f}s"                 

    def check_astrometry_data(self, header):
        return "CTYPE1" in header and "CTYPE2" in header

    def prompt_blind_solve(self):
        reply = QMessageBox.question(
            self, "Astrometry Data Missing",
            "No astrometry data found in the image. Would you like to perform a blind solve?",
            QMessageBox.Yes | QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            self.perform_blind_solve()

    def perform_blind_solve(self):
        # Load or prompt for API key
        api_key = load_api_key()
        if not api_key:
            api_key, ok = QInputDialog.getText(self, "Enter API Key", "Please enter your Astrometry.net API key:")
            if ok and api_key:
                save_api_key(api_key)
            else:
                QMessageBox.warning(self, "API Key Required", "Blind solve cannot proceed without an API key.")
                return

        try:
            self.status_label.setText("Status: Logging in to Astrometry.net...")
            QApplication.processEvents()

            # Step 1: Login to Astrometry.net
            session_key = self.login_to_astrometry(api_key)

            self.status_label.setText("Status: Uploading image to Astrometry.net...")
            QApplication.processEvents()
            
            # Step 2: Upload the image and get submission ID
            subid = self.upload_image_to_astrometry(self.image_path, session_key)

            self.status_label.setText("Status: Waiting for job ID...")
            QApplication.processEvents()
            
            # Step 3: Poll for the job ID until it's available
            job_id = self.poll_submission_status(subid)
            if not job_id:
                raise TimeoutError("Failed to retrieve job ID from Astrometry.net after multiple attempts.")
            
            self.status_label.setText("Status: Job ID found, processing image...")
            QApplication.processEvents()

            # Step 4a: Poll for the calibration data, ensuring RA/Dec are available
            calibration_data = self.poll_calibration_data(job_id)
            if not calibration_data:
                raise TimeoutError("Calibration data did not complete in the expected timeframe.")
            
            # Set pixscale and other necessary attributes from calibration data
            self.pixscale = calibration_data.get('pixscale')

            self.status_label.setText("Status: Calibration complete, downloading WCS file...")
            QApplication.processEvents()

            # Step 4b: Download the WCS FITS file for complete calibration data
            wcs_header = self.retrieve_and_apply_wcs(job_id)
            if not wcs_header:
                raise TimeoutError("Failed to retrieve WCS FITS file from Astrometry.net.")

            self.status_label.setText("Status: Applying astrometric solution to the image...")
            QApplication.processEvents()

            # Apply calibration data to the WCS
            self.apply_wcs_header(wcs_header)
            self.status_label.setText("Status: Blind Solve Complete.")
            QMessageBox.information(self, "Blind Solve Complete", "Astrometric solution applied successfully.")
        except Exception as e:
            self.status_label.setText("Status: Blind Solve Failed.")
            QMessageBox.critical(self, "Blind Solve Failed", f"An error occurred: {str(e)}")


    def retrieve_and_apply_wcs(self, job_id):
        """Download the wcs.fits file from Astrometry.net, extract WCS header data, and apply it."""
        try:
            wcs_url = f"https://nova.astrometry.net/wcs_file/{job_id}"
            wcs_filepath = "wcs.fits"
            max_retries = 10
            delay = 10  # seconds
            
            for attempt in range(max_retries):
                # Attempt to download the file
                response = requests.get(wcs_url, stream=True)
                response.raise_for_status()

                # Save the WCS file locally
                with open(wcs_filepath, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)

                # Check if the downloaded file is a valid FITS file
                try:
                    with fits.open(wcs_filepath, ignore_missing_simple=True, ignore_missing_end=True) as hdul:
                        # If it opens correctly, return the header
                        wcs_header = hdul[0].header
                        print("WCS header successfully retrieved.")
                        self.wcs = WCS(wcs_header)
                        return wcs_header
                except Exception as e:
                    print(f"Attempt {attempt + 1}: Failed to process WCS file - possibly HTML instead of FITS. Retrying in {delay} seconds...")
                    print(f"Error: {e}")
                    time.sleep(delay)  # Wait and retry
            
            print("Failed to download a valid WCS FITS file after multiple attempts.")
            return None

        except requests.exceptions.RequestException as e:
            print(f"Error downloading WCS file: {e}")
        except Exception as e:
            print(f"Error processing WCS file: {e}")
            
        return None



    def apply_wcs_header(self, wcs_header):
        """Apply WCS header to create a WCS object."""
        self.wcs = WCS(wcs_header)  # Initialize WCS with header directly
        print("WCS applied successfully from header data.")


    def calculate_pixel_from_ra_dec(self, ra, dec):
        """Convert RA/Dec to pixel coordinates using the WCS data."""
        if not hasattr(self, 'wcs'):
            print("WCS not initialized.")
            return None, None

        # Convert RA and Dec to pixel coordinates using the WCS object
        sky_coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
        x, y = self.wcs.world_to_pixel(sky_coord)
        
        return int(x), int(y)

    def login_to_astrometry(self, api_key):
        try:
            response = requests.post(
                ASTROMETRY_API_URL + "login",
                data={'request-json': json.dumps({"apikey": api_key})}
            )
            response_data = response.json()
            if response_data.get("status") == "success":
                return response_data["session"]
            else:
                raise ValueError("Login failed: " + response_data.get("error", "Unknown error"))
        except Exception as e:
            raise Exception("Login to Astrometry.net failed: " + str(e))

    def upload_image_to_astrometry(self, image_path, session_key):
        try:
            with open(image_path, 'rb') as image_file:
                files = {'file': image_file}
                data = {
                    'request-json': json.dumps({
                        "publicly_visible": "y",
                        "allow_modifications": "d",
                        "session": session_key,
                        "allow_commercial_use": "d"
                    })
                }
                response = requests.post(ASTROMETRY_API_URL + "upload", files=files, data=data)
                response_data = response.json()
                if response_data.get("status") == "success":
                    return response_data["subid"]
                else:
                    raise ValueError("Image upload failed: " + response_data.get("error", "Unknown error"))
        except Exception as e:
            raise Exception("Image upload to Astrometry.net failed: " + str(e))


    def poll_submission_status(self, subid):
        """Poll Astrometry.net to retrieve the job ID once the submission is processed."""
        max_retries = 90  # Adjust as necessary
        retries = 0
        while retries < max_retries:
            try:
                response = requests.get(ASTROMETRY_API_URL + f"submissions/{subid}")
                response_data = response.json()
                jobs = response_data.get("jobs", [])
                if jobs and jobs[0] is not None:
                    return jobs[0]
                else:
                    print(f"Polling attempt {retries + 1}: Job not ready yet.")
            except Exception as e:
                print(f"Error while polling submission status: {e}")
            
            retries += 1
            time.sleep(10)  # Wait 10 seconds between retries
        
        return None

    def poll_calibration_data(self, job_id):
        """Poll Astrometry.net to retrieve the calibration data once it's available."""
        max_retries = 90  # Retry for up to 15 minutes (90 * 10 seconds)
        retries = 0
        while retries < max_retries:
            try:
                response = requests.get(ASTROMETRY_API_URL + f"jobs/{job_id}/calibration/")
                response_data = response.json()
                if response_data and 'ra' in response_data and 'dec' in response_data:
                    print("Calibration data retrieved:", response_data)
                    return response_data  # Calibration data is complete
                else:
                    print(f"Calibration data not available yet (Attempt {retries + 1})")
            except Exception as e:
                print(f"Error retrieving calibration data: {e}")

            retries += 1
            time.sleep(10)  # Wait 10 seconds between retries

        return None


    #If originally a fits file update the header
    def update_fits_with_wcs(self, filepath, calibration_data):
        if not filepath.lower().endswith(('.fits', '.fit')):
            print("File is not a FITS file. Skipping WCS header update.")
            return

        print("Updating image with calibration data:", calibration_data)
        with fits.open(filepath, mode='update') as hdul:
            header = hdul[0].header
            header['CTYPE1'] = 'RA---TAN'
            header['CTYPE2'] = 'DEC--TAN'
            header['CRVAL1'] = calibration_data['ra']
            header['CRVAL2'] = calibration_data['dec']
            header['CRPIX1'] = hdul[0].data.shape[1] / 2
            header['CRPIX2'] = hdul[0].data.shape[0] / 2
            scale = calibration_data['pixscale'] / 3600
            orientation = np.radians(calibration_data['orientation'])
            header['CD1_1'] = -scale * np.cos(orientation)
            header['CD1_2'] = scale * np.sin(orientation)
            header['CD2_1'] = -scale * np.sin(orientation)
            header['CD2_2'] = -scale * np.cos(orientation)
            header['RADECSYS'] = 'ICRS'

    def on_mini_preview_press(self, event):
        # Set dragging flag and scroll the main preview to the position in the mini preview.
        self.dragging = True
        self.scroll_main_preview_to_mini_position(event)

    def on_mini_preview_drag(self, event):
        # Scroll to the new position while dragging in the mini preview.
        if self.dragging:
            self.scroll_main_preview_to_mini_position(event)

    def on_mini_preview_release(self, event):
        # Stop dragging
        self.dragging = False

    def scroll_main_preview_to_mini_position(self, event):
        """Scrolls the main preview to the corresponding position based on the mini preview click."""
        if self.main_image:
            # Get the click position in the mini preview
            click_x = event.pos().x()
            click_y = event.pos().y()
            
            # Calculate scale factors based on the difference in dimensions between main image and mini preview
            scale_factor_x = self.main_scene.sceneRect().width() / self.mini_preview.width()
            scale_factor_y = self.main_scene.sceneRect().height() / self.mini_preview.height()
            
            # Scale the click position to the main preview coordinates
            scaled_x = click_x * scale_factor_x
            scaled_y = click_y * scale_factor_y
            
            # Center the main preview on the calculated position
            self.main_preview.centerOn(scaled_x, scaled_y)
            
            # Update the green box after scrolling
            self.main_preview.update_mini_preview()

    def update_green_box(self):
        if self.main_image:
            factor_x = self.mini_preview.width() / self.main_image.width()
            factor_y = self.mini_preview.height() / self.main_image.height()
            
            # Get the current view rectangle in the main preview (in scene coordinates)
            view_rect = self.main_preview.mapToScene(self.main_preview.viewport().rect()).boundingRect()
            
            # Calculate the green box rectangle, shifted upward by half its height to center it
            green_box_rect = QRectF(
                view_rect.x() * factor_x,
                view_rect.y() * factor_y,
                view_rect.width() * factor_x,
                view_rect.height() * factor_y
            )
            
            # Scale the main image for the mini preview and draw the green box on it
            pixmap = self.main_image.scaled(self.mini_preview.size(), Qt.KeepAspectRatio, Qt.SmoothTransformation)
            painter = QPainter(pixmap)
            pen = QPen(QColor(0, 255, 0), 2)
            painter.setPen(pen)
            painter.drawRect(green_box_rect)
            painter.end()
            self.mini_preview.setPixmap(pixmap)



    def wheel_zoom(self, event):
        if event.angleDelta().y() > 0:
            self.zoom_in()
        else:
            self.zoom_out()

    def zoom_in(self):
        self.zoom_level *= 1.2
        self.main_preview.setTransform(QTransform().scale(self.zoom_level, self.zoom_level))
        self.update_green_box()
        
    def zoom_out(self):
        self.zoom_level /= 1.2
        self.main_preview.setTransform(QTransform().scale(self.zoom_level, self.zoom_level))
        self.update_green_box()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self.update_green_box()


    def update_circle_data(self):
        """Updates the status based on the circle's center and radius."""
        
        if self.circle_center and self.circle_radius > 0:
            if self.pixscale is None:
                print("Warning: Pixscale is None. Cannot calculate radius in arcminutes.")
                self.status_label.setText("No pixscale available for radius calculation.")
                return

            # Convert circle center to RA/Dec and radius to arcminutes
            ra, dec = self.calculate_ra_dec_from_pixel(self.circle_center.x(), self.circle_center.y())
            radius_arcmin = self.circle_radius * self.pixscale / 60.0  # Convert to arcminutes
            
            self.status_label.setText(
                f"Circle set at center RA={ra:.6f}, Dec={dec:.6f}, radius={radius_arcmin:.2f} arcmin"
            )
        else:
            self.status_label.setText("No search area defined.")

    def update_circle(self):
        """Draw the search circle on the main preview."""
        if self.circle_center and self.circle_radius > 0 and self.main_image:
            # Clear the main scene and add the main image back
            self.main_scene.clear()
            self.main_scene.addPixmap(self.main_image)

            # Draw the search circle
            pen_circle = QPen(QColor(255, 0, 0), 2)
            self.main_scene.addEllipse(
                int(self.circle_center.x() - self.circle_radius),
                int(self.circle_center.y() - self.circle_radius),
                int(self.circle_radius * 2),
                int(self.circle_radius * 2),
                pen_circle
            )

    def get_defined_radius(self):
        """Calculate radius in degrees for the defined region (circle radius)."""
        if self.circle_radius <= 0:
            return 0
        return float((self.circle_radius * self.pixscale) / 3600.0)


    def query_simbad(self, radius_deg, max_results=None):
        """Query Simbad based on the defined search circle using a single ADQL query, with filtering by selected types."""
            # If max_results is not provided, use the value from settings
        max_results = max_results if max_results is not None else self.max_results
        # Check if the circle center and radius are defined
        if not self.circle_center or self.circle_radius <= 0:
            QMessageBox.warning(self, "No Search Area", "Please define a search circle by Shift-clicking and dragging.")
            return

        # Calculate RA, Dec, and radius in degrees from pixel coordinates
        ra_center, dec_center = self.calculate_ra_dec_from_pixel(self.circle_center.x(), self.circle_center.y())
        if ra_center is None or dec_center is None:
            QMessageBox.warning(self, "Invalid Coordinates", "Could not determine the RA/Dec of the circle center.")
            return

        # Convert radius from arcseconds to degrees
        radius_deg = radius_deg

        # Get selected types from the tree widget
        selected_types = self.get_selected_object_types()
        if not selected_types:
            QMessageBox.warning(self, "No Object Types Selected", "Please select at least one object type.")
            return

        # Build ADQL query
        query = f"""
            SELECT TOP {max_results} ra, dec, main_id, rvz_redshift, otype, galdim_majaxis
            FROM basic
            WHERE CONTAINS(POINT('ICRS', basic.ra, basic.dec), CIRCLE('ICRS', {ra_center}, {dec_center}, {radius_deg})) = 1
        """

        try:
            # Execute the query using Simbad's TAP service
            result = Simbad.query_tap(query)

            # Clear previous results in the tree
            self.results_tree.clear()
            query_results = []

            if result is None or len(result) == 0:
                QMessageBox.information(self, "No Results", "No objects found in the specified area.")
                return

            # Process and display results, filtering by selected types
            for row in result:
                short_type = row["otype"]
                if short_type not in selected_types:
                    continue  # Skip items not in selected types

                # Retrieve other data fields
                ra = row["ra"]
                dec = row["dec"]
                main_id = row["main_id"]
                redshift = row["rvz_redshift"] if row["rvz_redshift"] is not None else "--"
                diameter = row.get("galdim_majaxis", "N/A")
                comoving_distance = calculate_comoving_distance(float(redshift)) if redshift != "--" else "N/A"

                # Map short type to long type
                long_type = otype_long_name_lookup.get(short_type, short_type)

                # Add to TreeWidget
                item = QTreeWidgetItem([
                    f"{ra:.6f}", f"{dec:.6f}", main_id, str(diameter), short_type, long_type, str(redshift), str(comoving_distance)
                ])
                self.results_tree.addTopLevelItem(item)

                # Append full details as a dictionary to query_results
                query_results.append({
                    'ra': ra,
                    'dec': dec,
                    'name': main_id,
                    'diameter': diameter,
                    'short_type': short_type,
                    'long_type': long_type,
                    'redshift': redshift,
                    'comoving_distance': comoving_distance,
                    'source' : "Simbad"
                })

            # Set query results in the CustomGraphicsView for display
            self.main_preview.set_query_results(query_results)
            self.query_results = query_results  # Keep a reference to results in MainWindow

        except Exception as e:
            QMessageBox.critical(self, "Query Failed", f"Failed to query Simbad: {str(e)}")

    def perform_deep_vizier_search(self):
        """Perform a Vizier catalog search and parse results based on catalog-specific fields, with duplicate handling."""
        if not self.circle_center or self.circle_radius <= 0:
            QMessageBox.warning(self, "No Search Area", "Please define a search circle by Shift-clicking and dragging.")
            return

        # Convert the center coordinates to RA/Dec
        ra_center, dec_center = self.calculate_ra_dec_from_pixel(self.circle_center.x(), self.circle_center.y())
        if ra_center is None or dec_center is None:
            QMessageBox.warning(self, "Invalid Coordinates", "Could not determine the RA/Dec of the circle center.")
            return

        # Convert radius from arcseconds to arcminutes
        radius_arcmin = float((self.circle_radius * self.pixscale) / 60.0)

        # List of Vizier catalogs
        catalog_ids = ["II/246", "I/350/gaiaedr3", "V/147/sdss12", "I/322A", "V/154"]

        coord = SkyCoord(ra_center, dec_center, unit="deg")
        all_results = []  # Collect all results for display in the main preview
        unique_entries = {}  # Dictionary to track unique entries by (RA, Dec) tuple

        try:
            for catalog_id in catalog_ids:
                # Query each catalog
                result = Vizier.query_region(coord, radius=radius_arcmin * u.arcmin, catalog=catalog_id)
                if result:
                    catalog_results = result[0]
                    for row in catalog_results:
                        # Map data to the columns in your tree view structure

                        # RA and Dec
                        ra = str(row.get("RAJ2000", row.get("RA_ICRS", "")))
                        dec = str(row.get("DEJ2000", row.get("DE_ICRS", "")))
                        if not ra or not dec:
                            
                            continue  # Skip this entry if RA or Dec is empty

                        # Create a unique key based on RA and Dec to track duplicates
                        unique_key = (ra, dec)

                        # Name (different columns based on catalog)
                        name = str(
                            row.get("_2MASS", "")
                            or row.get("Source", "")
                            or row.get("SDSS12", "")
                            or row.get("UCAC4", "")
                            or row.get("SDSS16", "")
                        )

                        # Diameter - store catalog ID as the diameter field to help with tracking
                        diameter = catalog_id

                        # Type (e.g., otype)
                        type_short = str(row.get("otype", "N/A"))

                        # Long Type (e.g., SpType)
                        long_type = str(row.get("SpType", "N/A"))

                        # Redshift or Parallax (zph for redshift or Plx for parallax)
                        redshift = row.get("zph", row.get("Plx", ""))
                        if redshift:
                            if "Plx" in row.colnames:
                                redshift = f"{redshift} (Parallax in mas)"
                                # Calculate the distance in light-years from parallax
                                try:
                                    parallax_value = float(row["Plx"])
                                    comoving_distance = f"{1000 / parallax_value * 3.2615637769:.2f} Ly"
                                except (ValueError, ZeroDivisionError):
                                    comoving_distance = "N/A"  # Handle invalid parallax values
                            else:
                                redshift = str(redshift)
                                # Calculate comoving distance for redshift if it's from zph
                                if "zph" in row.colnames and isinstance(row["zph"], (float, int)):
                                    comoving_distance = str(calculate_comoving_distance(float(row["zph"])))
                        else:
                            redshift = "N/A"
                            comoving_distance = "N/A"

                        # Handle duplicates: prioritize V/147/sdss12 over V/154 and only add unique entries
                        if unique_key not in unique_entries:
                            unique_entries[unique_key] = {
                                'ra': ra,
                                'dec': dec,
                                'name': name,
                                'diameter': diameter,
                                'short_type': type_short,
                                'long_type': long_type,
                                'redshift': redshift,
                                'comoving_distance': comoving_distance,
                                'source' : "Vizier"
                            }
                        else:
                            # Check if we should replace the existing entry
                            existing_entry = unique_entries[unique_key]
                            if (existing_entry['diameter'] == "V/154" and diameter == "V/147/sdss12"):
                                unique_entries[unique_key] = {
                                    'ra': ra,
                                    'dec': dec,
                                    'name': name,
                                    'diameter': diameter,
                                    'short_type': type_short,
                                    'long_type': long_type,
                                    'redshift': redshift,
                                    'comoving_distance': comoving_distance,
                                    'source' : "Vizier"
                                }

            # Convert unique entries to the main preview display
            for entry in unique_entries.values():
                item = QTreeWidgetItem([
                    entry['ra'], entry['dec'], entry['name'], entry['diameter'], entry['short_type'], entry['long_type'],
                    entry['redshift'], entry['comoving_distance']
                ])
                self.results_tree.addTopLevelItem(item)
                all_results.append(entry)

            # Update the main preview with the query results
            self.main_preview.set_query_results(all_results)
            self.query_results = all_results  # Keep a reference to results in MainWindow
            
        except Exception as e:
            QMessageBox.critical(self, "Vizier Search Failed", f"Failed to query Vizier: {str(e)}")




    def toggle_show_names(self, state):
        """Toggle showing/hiding names on the main image."""
        self.show_names = state == Qt.Checked
        self.main_preview.draw_query_results()  # Redraw with or without names

    def clear_results(self):
        """Clear the search results and remove markers from the main image."""
        self.results_tree.clear()
        self.main_preview.clear_query_results()
        self.status_label.setText("Results cleared.")

    def open_settings_dialog(self):
        """Open settings dialog to adjust max results."""
        dialog = QDialog(self)
        dialog.setWindowTitle("Settings")
        
        layout = QFormLayout(dialog)
        
        # Max Results setting
        max_results_spinbox = QSpinBox()
        max_results_spinbox.setRange(1, 25000)
        max_results_spinbox.setValue(self.max_results)
        layout.addRow("Max Results:", max_results_spinbox)
        
        # OK and Cancel buttons
        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(lambda: self.update_max_results(max_results_spinbox.value(), dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        
        dialog.setLayout(layout)
        dialog.exec_()

    def update_max_results(self, value, dialog):
        """Update the max_results variable and close the dialog."""
        self.max_results = value
        dialog.accept()            


def extract_wcs_data(file_path):
    try:
        # Open the FITS file with minimal validation to ignore potential errors in non-essential parts
        with fits.open(file_path, ignore_missing_simple=True, ignore_missing_end=True) as hdul:
            header = hdul[0].header

            # Extract essential WCS parameters
            wcs_params = {}
            keys_to_extract = [
                'WCSAXES', 'CTYPE1', 'CTYPE2', 'EQUINOX', 'LONPOLE', 'LATPOLE',
                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2', 'CUNIT1', 'CUNIT2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'A_ORDER', 'A_0_0', 'A_0_1', 
                'A_0_2', 'A_1_0', 'A_1_1', 'A_2_0', 'B_ORDER', 'B_0_0', 'B_0_1', 
                'B_0_2', 'B_1_0', 'B_1_1', 'B_2_0', 'AP_ORDER', 'AP_0_0', 'AP_0_1', 
                'AP_0_2', 'AP_1_0', 'AP_1_1', 'AP_2_0', 'BP_ORDER', 'BP_0_0', 
                'BP_0_1', 'BP_0_2', 'BP_1_0', 'BP_1_1', 'BP_2_0'
            ]
            for key in keys_to_extract:
                if key in header:
                    wcs_params[key] = header[key]

            # Manually create a minimal header with WCS information
            wcs_header = fits.Header()
            for key, value in wcs_params.items():
                wcs_header[key] = value

            # Initialize WCS with this custom header
            wcs = WCS(wcs_header)
            print("WCS successfully initialized with minimal header.")
            return wcs

    except Exception as e:
        print(f"Error processing WCS file: {e}")
        return None

# Function to calculate comoving radial distance (in Gly)
def calculate_comoving_distance(z):
    z = abs(z)
    # Initialize variables
    WR = 4.165E-5 / ((H0 / 100) ** 2)  # Omega radiation
    WK = 1 - WM - WV - WR  # Omega curvature
    az = 1.0 / (1 + z)
    n = 1000  # number of points in integration

    # Comoving radial distance
    DCMR = 0.0
    for i in range(n):
        a = az + (1 - az) * (i + 0.5) / n
        adot = sqrt(WK + (WM / a) + (WR / (a ** 2)) + (WV * a ** 2))
        DCMR += 1 / (a * adot)
    
    DCMR = (1 - az) * DCMR / n
    DCMR_Gly = (c / H0) * DCMR * Mpc_to_Gly

    return round(DCMR_Gly, 3)  # Round to three decimal places for display

app = QApplication(sys.argv)
window = MainWindow()
window.show()
sys.exit(app.exec_())
