import sys
import os
import numpy as np
import tifffile as tiff
from PIL import Image
from astropy.io import fits
from PyQt5.QtWidgets import (QApplication, QWidget, QTabWidget, QVBoxLayout, QLabel, QPushButton, 
                             QSlider, QCheckBox, QFileDialog, QHBoxLayout, QSpacerItem, 
                             QSizePolicy, QScrollArea, QInputDialog, QComboBox, QRadioButton, QGridLayout)
from PyQt5.QtCore import Qt, QPoint, QThread, pyqtSignal, QObject, QTimer, QCoreApplication
from PyQt5.QtGui import QPixmap, QImage, QMovie
import matplotlib.colors as colors  # For color conversions
import cv2


class AstroEditingSuite(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        layout = QVBoxLayout()

        # Create tab widget
        self.tabs = QTabWidget()

        # Add individual tabs for each tool
        self.tabs.addTab(StatisticalStretchTab(), "Statistical Stretch")
        self.tabs.addTab(NBtoRGBstarsTab(), "NB to RGB Stars")  # Placeholder        
        self.tabs.addTab(StarStretchTab(), "Star Stretch")  # Placeholder
        self.tabs.addTab(HaloBGonTab(), "Halo-B-Gon")  # Placeholder
        self.tabs.addTab(ContinuumSubtractTab(), "Continuum Subtraction")

        # Add the tab widget to the main layout
        layout.addWidget(self.tabs)

        # Set the layout for the main window
        self.setLayout(layout)
        self.setWindowTitle('Seti Astro\'s Editing Suite V1.0')


class StatisticalStretchTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.image = None
        self.filename = None
        self.zoom_factor = 1.0
        self.original_header = None

    def initUI(self):
        main_layout = QHBoxLayout()

        # Left column for controls
        left_widget = QWidget(self)
        left_layout = QVBoxLayout(left_widget)
        left_widget.setFixedWidth(400)  # You can adjust this width as needed

        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Select an image to stretch.
            2. Adjust the target median and optional settings.
            3. Preview the result.
            4. Save the stretched image in your desired format.
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # File selection button
        self.fileButton = QPushButton('Select Image', self)
        self.fileButton.clicked.connect(self.openFileDialog)
        left_layout.addWidget(self.fileButton)

        self.fileLabel = QLabel('No file selected', self)
        left_layout.addWidget(self.fileLabel)

        # Target median slider
        self.medianLabel = QLabel('Target Median: 0.25', self)
        self.medianSlider = QSlider(Qt.Horizontal)
        self.medianSlider.setMinimum(1)
        self.medianSlider.setMaximum(100)
        self.medianSlider.setValue(25)
        self.medianSlider.valueChanged.connect(self.updateMedianLabel)
        left_layout.addWidget(self.medianLabel)
        left_layout.addWidget(self.medianSlider)

        # Linked/Unlinked stretch checkbox
        self.linkedCheckBox = QCheckBox('Linked Stretch', self)
        self.linkedCheckBox.setChecked(True)
        left_layout.addWidget(self.linkedCheckBox)

        # Normalization checkbox
        self.normalizeCheckBox = QCheckBox('Normalize Image', self)
        left_layout.addWidget(self.normalizeCheckBox)

        # Curves adjustment checkbox
        self.curvesCheckBox = QCheckBox('Apply Curves Adjustment', self)
        self.curvesCheckBox.stateChanged.connect(self.toggleCurvesSlider)
        left_layout.addWidget(self.curvesCheckBox)

        # Curves Boost slider (initially hidden)
        self.curvesBoostLabel = QLabel('Curves Boost: 0.00', self)
        self.curvesBoostSlider = QSlider(Qt.Horizontal)
        self.curvesBoostSlider.setMinimum(0)
        self.curvesBoostSlider.setMaximum(50)
        self.curvesBoostSlider.setValue(0)
        self.curvesBoostSlider.valueChanged.connect(self.updateCurvesBoostLabel)
        self.curvesBoostLabel.hide()
        self.curvesBoostSlider.hide()

        left_layout.addWidget(self.curvesBoostLabel)
        left_layout.addWidget(self.curvesBoostSlider)

        # Progress indicator (spinner) label
        self.spinnerLabel = QLabel(self)
        self.spinnerLabel.setAlignment(Qt.AlignCenter)
        # Use the resource path function to access the GIF
        self.spinnerMovie = QMovie(resource_path("spinner.gif"))  # Updated path
        self.spinnerLabel.setMovie(self.spinnerMovie)
        self.spinnerLabel.hide()  # Hide spinner by default
        left_layout.addWidget(self.spinnerLabel)      

        # Preview button
        self.previewButton = QPushButton('Preview Stretch', self)
        self.previewButton.clicked.connect(self.previewStretch)
        left_layout.addWidget(self.previewButton)

        # Zoom buttons
        zoom_layout = QHBoxLayout()
        self.zoomInButton = QPushButton('Zoom In', self)
        self.zoomInButton.clicked.connect(self.zoom_in)
        zoom_layout.addWidget(self.zoomInButton)

        self.zoomOutButton = QPushButton('Zoom Out', self)
        self.zoomOutButton.clicked.connect(self.zoom_out)
        zoom_layout.addWidget(self.zoomOutButton)

        left_layout.addLayout(zoom_layout)

        # Save button
        self.saveButton = QPushButton('Save Stretched Image', self)
        self.saveButton.clicked.connect(self.saveImage)
        left_layout.addWidget(self.saveButton)

                        # Footer
        footer_label = QLabel("""
            Written by Franklin Marek<br>
            <a href='http://www.setiastro.com'>www.setiastro.com</a>
        """)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setOpenExternalLinks(True)
        footer_label.setStyleSheet("font-size: 10px;")
        left_layout.addWidget(footer_label)

        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Add the left widget to the main layout
        main_layout.addWidget(left_widget)

        # Right side for the preview inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        # QLabel for the image preview
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400, 400)

        main_layout.addWidget(self.scrollArea)

        self.setLayout(main_layout)
        self.zoom_factor = 0.25
        self.scrollArea.viewport().setMouseTracking(True)
        self.scrollArea.viewport().installEventFilter(self)
        self.dragging = False
        self.last_pos = QPoint()

    def eventFilter(self, source, event):
        if event.type() == event.MouseButtonPress and event.button() == Qt.LeftButton:
            self.dragging = True
            self.last_pos = event.pos()
        elif event.type() == event.MouseButtonRelease and event.button() == Qt.LeftButton:
            self.dragging = False
        elif event.type() == event.MouseMove and self.dragging:
            delta = event.pos() - self.last_pos
            self.scrollArea.horizontalScrollBar().setValue(self.scrollArea.horizontalScrollBar().value() - delta.x())
            self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().value() - delta.y())
            self.last_pos = event.pos()

        return super().eventFilter(source, event)


    def openFileDialog(self):
        self.filename, _ = QFileDialog.getOpenFileName(self, 'Open Image File', '', 'Images (*.png *.tiff *.tif *.fits *.fit);;All Files (*)')
        if self.filename:
            self.fileLabel.setText(self.filename)
            # Unpack all four values returned by load_image
            self.image, self.original_header, self.bit_depth, self.is_mono = load_image(self.filename)


    def updateMedianLabel(self, value):
        self.medianLabel.setText(f'Target Median: {value / 100:.2f}')

    def updateCurvesBoostLabel(self, value):
        self.curvesBoostLabel.setText(f'Curves Boost: {value / 100:.2f}')

    def toggleCurvesSlider(self, state):
        if state == Qt.Checked:
            self.curvesBoostLabel.show()
            self.curvesBoostSlider.show()
        else:
            self.curvesBoostLabel.hide()
            self.curvesBoostSlider.hide()

    def previewStretch(self):
        if self.image is not None:
            # Show spinner before starting processing
            self.showSpinner()

            # Start background processing
            self.processing_thread = StatisticalStretchProcessingThread(self.image,
                                                                        self.medianSlider.value(),
                                                                        self.linkedCheckBox.isChecked(),
                                                                        self.normalizeCheckBox.isChecked(),
                                                                        self.curvesCheckBox.isChecked(),
                                                                        self.curvesBoostSlider.value() / 100.0)
            self.processing_thread.preview_generated.connect(self.update_preview)
            self.processing_thread.start()

    def update_preview(self, stretched_image):
        # Save the stretched image for later use in zoom functions
        self.stretched_image = stretched_image

        # Update the preview once the processing thread emits the result
        img = (stretched_image * 255).astype(np.uint8)
        h, w = img.shape[:2]

        if img.ndim == 3:
            bytes_per_line = 3 * w
            q_image = QImage(img.tobytes(), w, h, bytes_per_line, QImage.Format_RGB888)
        else:
            bytes_per_line = w
            q_image = QImage(img.tobytes(), w, h, bytes_per_line, QImage.Format_Grayscale8)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())

        # Hide the spinner after processing is done
        self.hideSpinner()


    def showSpinner(self):
        self.spinnerLabel.show()
        self.spinnerMovie.start()

    def hideSpinner(self):
        self.spinnerLabel.hide()
        self.spinnerMovie.stop()    

    def zoom_in(self):
        if hasattr(self, 'stretched_image'):
            self.zoom_factor *= 1.2
            self.update_preview(self.stretched_image)  # Pass the latest stretched image

    def zoom_out(self):
        if hasattr(self, 'stretched_image'):
            self.zoom_factor /= 1.2
            self.update_preview(self.stretched_image)  # Pass the latest stretched image



    def saveImage(self):
        if self.image is not None:
            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename)
            default_save_name = os.path.splitext(base_name)[0] + '_stretched.tif'
            original_dir = os.path.dirname(self.filename)

            # Open the save file dialog
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name), 
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Call save_image with the necessary parameters
                        save_image(self.image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # For non-TIFF/FITS formats, save directly without bit depth selection
                    save_image(self.image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')


# Thread for Stat Stretch background processing
class StatisticalStretchProcessingThread(QThread):
    preview_generated = pyqtSignal(np.ndarray)  # Signal to send the generated preview image back to the main thread

    def __init__(self, image, target_median, linked, normalize, apply_curves, curves_boost):
        super().__init__()
        self.image = image
        self.target_median = target_median / 100.0  # Ensure proper scaling
        self.linked = linked
        self.normalize = normalize
        self.apply_curves = apply_curves
        self.curves_boost = curves_boost

    def run(self):
        # Perform the image stretching in the background
        if self.image.ndim == 2:  # Mono image
            stretched_image = stretch_mono_image(self.image, self.target_median, self.normalize, self.apply_curves, self.curves_boost)
        else:  # Color image
            stretched_image = stretch_color_image(self.image, self.target_median, self.linked, self.normalize, self.apply_curves, self.curves_boost)

        # Emit the result once done
        self.preview_generated.emit(stretched_image)

# Thread for star stretch background processing
class ProcessingThread(QThread):
    preview_generated = pyqtSignal(np.ndarray)

    def __init__(self, image, stretch_factor, sat_amount, scnr_enabled):
        super().__init__()
        self.image = image
        self.stretch_factor = stretch_factor
        self.sat_amount = sat_amount
        self.scnr_enabled = scnr_enabled

    def run(self):
        stretched_image = self.applyPixelMath(self.image, self.stretch_factor)
        stretched_image = self.applyColorSaturation(stretched_image, self.sat_amount)
        if self.scnr_enabled:
            stretched_image = self.applySCNR(stretched_image)
        self.preview_generated.emit(stretched_image)

    def applyPixelMath(self, image_array, amount):
        expression = (3 ** amount * image_array) / ((3 ** amount - 1) * image_array + 1)
        return np.clip(expression, 0, 1)

    def applyColorSaturation(self, image_array, satAmount):
        saturationLevel = [
            [0.0, satAmount * 0.4],
            [0.5, satAmount * 0.7],
            [1.0, satAmount * 0.4]
        ]
        return self.adjust_saturation(image_array, saturationLevel)

    def adjust_saturation(self, image_array, saturation_level):
        hsv_image = np.array(Image.fromarray((image_array * 255).astype(np.uint8)).convert('HSV')) / 255.0
        hsv_image[..., 1] *= saturation_level[1][1]
        hsv_image[..., 1] = np.clip(hsv_image[..., 1], 0, 1)
        rgb_image = Image.fromarray((hsv_image * 255).astype(np.uint8), 'HSV').convert('RGB')
        return np.array(rgb_image) / 255.0

    def applySCNR(self, image_array):
        red_channel = image_array[..., 0]
        green_channel = image_array[..., 1]
        blue_channel = image_array[..., 2]

        # Apply green neutralization where green is higher than red and blue
        mask = green_channel > np.maximum(red_channel, blue_channel)
        green_channel[mask] = np.maximum(red_channel[mask], blue_channel[mask])

        # Recombine the channels
        image_array[..., 1] = green_channel
        return np.clip(image_array, 0, 1)


class StarStretchTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.image = None  # Store the selected image
        self.stretch_factor = 5.0
        self.sat_amount = 1.0
        self.is_mono = True
        self.remove_green = False
        self.filename = None  # Store the selected file path
        self.preview_image = None  # Store the preview result
        self.zoom_factor = 0.25  # Initialize zoom factor for preview scaling
        self.dragging = False
        self.last_pos = None
        self.processing_thread = None  # Thread for background processing
        self.original_header = None

    def initUI(self):
        main_layout = QHBoxLayout()

        # Left column for controls
        left_widget = QWidget(self)
        left_layout = QVBoxLayout(left_widget)
        left_widget.setFixedWidth(400)  # Fix the left column width

        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Select a stars-only image.
            2. Adjust the stretch and optional settings.
            3. Preview the result.
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # File selection button
        self.fileButton = QPushButton("Select Stars Only Image", self)
        self.fileButton.clicked.connect(self.selectImage)
        left_layout.addWidget(self.fileButton)

        self.fileLabel = QLabel('No file selected', self)
        left_layout.addWidget(self.fileLabel)

        # Stretch Amount slider with more precision
        self.stretchLabel = QLabel("Stretch Amount: 5.00", self)
        self.stretchSlider = QSlider(Qt.Horizontal)
        self.stretchSlider.setMinimum(0)
        self.stretchSlider.setMaximum(800)  # Allow two decimal places of precision
        self.stretchSlider.setValue(500)  # 500 corresponds to 5.00
        self.stretchSlider.valueChanged.connect(self.updateStretchLabel)
        left_layout.addWidget(self.stretchLabel)
        left_layout.addWidget(self.stretchSlider)

        # Color Boost Amount slider
        self.satLabel = QLabel("Color Boost: 1.00", self)
        self.satSlider = QSlider(Qt.Horizontal)
        self.satSlider.setMinimum(0)
        self.satSlider.setMaximum(200)
        self.satSlider.setValue(100)  # 100 corresponds to 1.0 boost
        self.satSlider.valueChanged.connect(self.updateSatLabel)
        left_layout.addWidget(self.satLabel)
        left_layout.addWidget(self.satSlider)

        # SCNR checkbox
        self.scnrCheckBox = QCheckBox("Remove Green via SCNR (Optional)", self)
        left_layout.addWidget(self.scnrCheckBox)

        # Refresh preview button
        self.refreshButton = QPushButton("Refresh Preview", self)
        self.refreshButton.clicked.connect(self.generatePreview)
        left_layout.addWidget(self.refreshButton)

        # Progress indicator (spinner) label
        self.spinnerLabel = QLabel(self)
        self.spinnerLabel.setAlignment(Qt.AlignCenter)
        # Use the resource path function to access the GIF
        self.spinnerMovie = QMovie(resource_path("spinner.gif"))  # Updated path
        self.spinnerLabel.setMovie(self.spinnerMovie)
        self.spinnerLabel.hide()  # Hide spinner by default
        left_layout.addWidget(self.spinnerLabel)

        # Zoom buttons for preview
        zoom_layout = QHBoxLayout()
        self.zoomInButton = QPushButton('Zoom In', self)
        self.zoomInButton.clicked.connect(self.zoom_in)
        zoom_layout.addWidget(self.zoomInButton)

        self.zoomOutButton = QPushButton('Zoom Out', self)
        self.zoomOutButton.clicked.connect(self.zoom_out)
        zoom_layout.addWidget(self.zoomOutButton)
        left_layout.addLayout(zoom_layout)

        # Save As button (replaces Execute button)
        self.saveAsButton = QPushButton("Save As", self)
        self.saveAsButton.clicked.connect(self.saveImage)
        left_layout.addWidget(self.saveAsButton)

                        # Footer
        footer_label = QLabel("""
            Written by Franklin Marek<br>
            <a href='http://www.setiastro.com'>www.setiastro.com</a>
        """)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setOpenExternalLinks(True)
        footer_label.setStyleSheet("font-size: 10px;")
        left_layout.addWidget(footer_label)

        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Add the left widget to the main layout
        main_layout.addWidget(left_widget)

        # Right side for the preview inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.viewport().installEventFilter(self)

        # QLabel for the image preview
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400, 400)

        main_layout.addWidget(self.scrollArea)

        self.setLayout(main_layout)

    def saveImage(self):
        if self.image is not None:
            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename)
            default_save_name = os.path.splitext(base_name)[0] + '_stretched.tif'
            original_dir = os.path.dirname(self.filename)

            # Open the save file dialog
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name), 
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Call save_image with the necessary parameters
                        save_image(self.image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # For non-TIFF/FITS formats, save directly without bit depth selection
                    save_image(self.image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')


    def selectImage(self):
        selected_file, _ = QFileDialog.getOpenFileName(self, "Select Stars Only Image", "", "Images (*.png *.tif *.fits *.fit)")
        if selected_file:
            try:
                self.image, self.original_header, _, self.is_mono = load_image(selected_file)  # Load image with header
                self.filename = selected_file  # Store the selected file path
                self.fileLabel.setText(os.path.basename(selected_file))
                self.generatePreview()

            except Exception as e:
                self.fileLabel.setText(f"Error: {str(e)}")
                print(f"Failed to load image: {e}")

    def updateStretchLabel(self, value):
        self.stretch_factor = value / 100.0  # Precision of two decimals
        self.stretchLabel.setText(f"Stretch Amount: {self.stretch_factor:.2f}")

    def updateSatLabel(self, value):
        self.sat_amount = value / 100.0
        self.satLabel.setText(f"Color Boost: {self.sat_amount:.2f}")

    def generatePreview(self):
        if self.image is not None and self.image.size > 0:
            # Show spinner before starting processing
            self.showSpinner()

            # Start background processing
            self.processing_thread = ProcessingThread(self.image, self.stretch_factor, self.sat_amount, self.scnrCheckBox.isChecked())
            self.processing_thread.preview_generated.connect(self.updatePreview)
            self.processing_thread.start()

    def updatePreview(self, stretched_image):
        # Update the preview once the processing thread emits the result
        preview_image = (stretched_image * 255).astype(np.uint8)
        h, w = preview_image.shape[:2]
        if preview_image.ndim == 3:
            q_image = QImage(preview_image.data, w, h, 3 * w, QImage.Format_RGB888)
        else:
            q_image = QImage(preview_image.data, w, h, w, QImage.Format_Grayscale8)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())

        # Hide the spinner after processing is done
        self.hideSpinner()

    def eventFilter(self, source, event):
        if event.type() == event.MouseButtonPress and event.button() == Qt.LeftButton:
            self.dragging = True
            self.last_pos = event.pos()
        elif event.type() == event.MouseButtonRelease and event.button() == Qt.LeftButton:
            self.dragging = False
        elif event.type() == event.MouseMove and self.dragging:
            delta = event.pos() - self.last_pos
            self.scrollArea.horizontalScrollBar().setValue(self.scrollArea.horizontalScrollBar().value() - delta.x())
            self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().value() - delta.y())
            self.last_pos = event.pos()

        return super().eventFilter(source, event)
    

    def showSpinner(self):
        self.spinnerLabel.show()
        self.spinnerMovie.start()

    def hideSpinner(self):
        self.spinnerLabel.hide()
        self.spinnerMovie.stop()

    def zoom_in(self):
        self.zoom_factor *= 1.2
        self.generatePreview()

    def zoom_out(self):
        self.zoom_factor /= 1.2
        self.generatePreview()

    def applyStretch(self):
        if self.image is not None and self.image.size > 0:
            print(f"Applying stretch: {self.stretch_factor}, Color Boost: {self.sat_amount:.2f}, SCNR: {self.scnrCheckBox.isChecked()}")
            self.generatePreview()

class NBtoRGBstarsTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.ha_image = None
        self.oiii_image = None
        self.sii_image = None
        self.osc_image = None
        self.combined_image = None
        self.is_mono = False
        self.filename = None  # Store the selected file path
        self.zoom_factor = 0.25
        self.dragging = False
        self.last_pos = QPoint()
        self.processing_thread = None
        self.original_header = None

    def initUI(self):
        main_layout = QHBoxLayout()

        # Left column for controls
        left_widget = QWidget(self)
        left_layout = QVBoxLayout(left_widget)
        left_widget.setFixedWidth(400)

        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Select Ha, OIII, and SII (optional) narrowband images, or an OSC stars-only image.
            2. Adjust the Ha to OIII Ratio if needed.
            3. Preview the combined result.
            4. Save the final composite image.
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # Ha, OIII, SII image selections
        self.haButton = QPushButton('Select Ha Image', self)
        self.haButton.clicked.connect(lambda: self.selectImage('Ha'))
        left_layout.addWidget(self.haButton)
        self.haLabel = QLabel('No Ha image selected', self)
        left_layout.addWidget(self.haLabel)

        self.oiiiButton = QPushButton('Select OIII Image', self)
        self.oiiiButton.clicked.connect(lambda: self.selectImage('OIII'))
        left_layout.addWidget(self.oiiiButton)
        self.oiiiLabel = QLabel('No OIII image selected', self)
        left_layout.addWidget(self.oiiiLabel)

        self.siiButton = QPushButton('Select SII Image (Optional)', self)
        self.siiButton.clicked.connect(lambda: self.selectImage('SII'))
        left_layout.addWidget(self.siiButton)
        self.siiLabel = QLabel('No SII image selected', self)
        left_layout.addWidget(self.siiLabel)

        self.oscButton = QPushButton('Select OSC Stars Image (Optional)', self)
        self.oscButton.clicked.connect(lambda: self.selectImage('OSC'))
        left_layout.addWidget(self.oscButton)
        self.oscLabel = QLabel('No OSC stars image selected', self)
        left_layout.addWidget(self.oscLabel)

        # Ha to OIII Ratio slider
        self.haToOiiRatioLabel, self.haToOiiRatioSlider = self.createRatioSlider("Ha to OIII Ratio", 30)
        left_layout.addWidget(self.haToOiiRatioLabel)
        left_layout.addWidget(self.haToOiiRatioSlider)

        # Star Stretch checkbox and sliders
        self.starStretchCheckBox = QCheckBox("Enable Star Stretch", self)
        self.starStretchCheckBox.setChecked(True)
        self.starStretchCheckBox.toggled.connect(self.toggleStarStretchControls)
        left_layout.addWidget(self.starStretchCheckBox)

        self.stretchSliderLabel, self.stretchSlider = self.createStretchSlider("Stretch Factor", 5.0)
        left_layout.addWidget(self.stretchSliderLabel)
        left_layout.addWidget(self.stretchSlider)

        # Progress indicator (spinner) label
        self.spinnerLabel = QLabel(self)
        self.spinnerLabel.setAlignment(Qt.AlignCenter)
        # Use the resource path function to access the GIF
        self.spinnerMovie = QMovie(resource_path("spinner.gif"))  # Updated path
        self.spinnerLabel.setMovie(self.spinnerMovie)
        self.spinnerLabel.hide()  # Hide spinner by default
        left_layout.addWidget(self.spinnerLabel)

        # Preview and Save buttons
        self.previewButton = QPushButton('Preview Combined Image', self)
        self.previewButton.clicked.connect(self.previewCombine)
        left_layout.addWidget(self.previewButton)

        # File label for displaying save status
        self.fileLabel = QLabel('', self)
        left_layout.addWidget(self.fileLabel)

        self.saveButton = QPushButton('Save Combined Image', self)
        self.saveButton.clicked.connect(self.saveImage)
        left_layout.addWidget(self.saveButton)

                        # Footer
        footer_label = QLabel("""
            Written by Franklin Marek<br>
            <a href='http://www.setiastro.com'>www.setiastro.com</a>
        """)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setOpenExternalLinks(True)
        footer_label.setStyleSheet("font-size: 10px;")
        left_layout.addWidget(footer_label)

        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))
        main_layout.addWidget(left_widget)

        # Right side for the preview inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400, 400)
        main_layout.addWidget(self.scrollArea)

        self.setLayout(main_layout)
        self.scrollArea.viewport().setMouseTracking(True)
        self.scrollArea.viewport().installEventFilter(self)

    def createRatioSlider(self, label_text, default_value):
        label = QLabel(f"{label_text}: {default_value / 100:.2f}", self)
        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(100)
        slider.setValue(default_value)
        slider.valueChanged.connect(lambda value: label.setText(f"{label_text}: {value / 100:.2f}"))
        return label, slider

    def createStretchSlider(self, label_text, default_value):
        label = QLabel(f"{label_text}: {default_value:.2f}", self)
        slider = QSlider(Qt.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(800)
        slider.setValue(int(default_value * 100))  # Scale to handle float values
        slider.valueChanged.connect(lambda value: label.setText(f"{label_text}: {value / 100:.2f}"))
        return label, slider

    def toggleStarStretchControls(self):
        enabled = self.starStretchCheckBox.isChecked()
        self.stretchSliderLabel.setVisible(enabled)
        self.stretchSlider.setVisible(enabled)

    def selectImage(self, image_type):
        selected_file, _ = QFileDialog.getOpenFileName(self, f"Select {image_type} Image", "", "Images (*.png *.tif *.fits *.fit)")
        if selected_file:
            try:
                if image_type == 'Ha':
                    self.ha_image, self.original_header, _, _ = load_image(selected_file)  # Store header
                    self.filename = selected_file
                    self.haLabel.setText(f"{os.path.basename(selected_file)} selected")
                elif image_type == 'OIII':
                    self.oiii_image, self.original_header, _, _ = load_image(selected_file)
                    self.oiiiLabel.setText(f"{os.path.basename(selected_file)} selected")
                elif image_type == 'SII':
                    self.sii_image, self.original_header, _, _ = load_image(selected_file)
                    self.siiLabel.setText(f"{os.path.basename(selected_file)} selected")
                elif image_type == 'OSC':
                    self.osc_image, self.original_header, _, _ = load_image(selected_file)
                    self.oscLabel.setText(f"{os.path.basename(selected_file)} selected")
            except Exception as e:
                print(f"Failed to load {image_type} image: {e}")


    def previewCombine(self):
        ha_to_oii_ratio = self.haToOiiRatioSlider.value() / 100.0
        enable_star_stretch = self.starStretchCheckBox.isChecked()
        stretch_factor = self.stretchSlider.value() / 100.0

        # Show spinner before starting processing
        self.showSpinner()

        # Start background processing
        self.processing_thread = NBtoRGBProcessingThread(
            self.ha_image, self.oiii_image, self.sii_image, self.osc_image,
            ha_to_oii_ratio=ha_to_oii_ratio, enable_star_stretch=enable_star_stretch, stretch_factor=stretch_factor
        )
        self.processing_thread.preview_generated.connect(self.updatePreview)
        self.processing_thread.start()

    def updatePreview(self, combined_image):
        # Set the combined image for saving
        self.combined_image = combined_image

        # Convert the image to display format
        preview_image = (combined_image * 255).astype(np.uint8)
        h, w = preview_image.shape[:2]
        q_image = QImage(preview_image.data, w, h, 3 * w, QImage.Format_RGB888)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())

        # Hide the spinner after processing is done
        self.hideSpinner()

    def showSpinner(self):
        self.spinnerLabel.show()
        self.spinnerMovie.start()

    def hideSpinner(self):
        self.spinnerLabel.hide()
        self.spinnerMovie.stop()

    def saveImage(self):
        if self.combined_image is not None:
            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename) if self.filename else "output"
            default_save_name = 'NBtoRGBstars.tif'
            original_dir = os.path.dirname(self.filename) if self.filename else ""

            # Open the save file dialog
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name), 
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Call save_image with the necessary parameters
                        save_image(self.combined_image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # For non-TIFF/FITS formats, save directly without bit depth selection
                    save_image(self.combined_image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')
        else:
            self.fileLabel.setText("No combined image to save.")




    # Add event filter for mouse dragging in preview area
    def eventFilter(self, source, event):
        if event.type() == event.MouseButtonPress and event.button() == Qt.LeftButton:
            self.dragging = True
            self.last_pos = event.pos()
        elif event.type() == event.MouseButtonRelease and event.button() == Qt.LeftButton:
            self.dragging = False
        elif event.type() == event.MouseMove and self.dragging:
            delta = event.pos() - self.last_pos
            self.scrollArea.horizontalScrollBar().setValue(self.scrollArea.horizontalScrollBar().value() - delta.x())
            self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().value() - delta.y())
            self.last_pos = event.pos()

        return super().eventFilter(source, event)

class NBtoRGBProcessingThread(QThread):
    preview_generated = pyqtSignal(np.ndarray)

    def __init__(self, ha_image, oiii_image, sii_image=None, osc_image=None, ha_to_oii_ratio=0.3, enable_star_stretch=True, stretch_factor=5.0):
        super().__init__()
        self.ha_image = ha_image
        self.oiii_image = oiii_image
        self.sii_image = sii_image
        self.osc_image = osc_image
        self.ha_to_oii_ratio = ha_to_oii_ratio
        self.enable_star_stretch = enable_star_stretch
        self.stretch_factor = stretch_factor

    def run(self):
        # Check if an OSC image is provided and handle as RGB channels
        if self.osc_image is not None:
            r_channel = self.osc_image[..., 0]
            g_channel = self.osc_image[..., 1]
            b_channel = self.osc_image[..., 2]
            
            # Handle cases where narrowband images are missing
            # If Ha is None, use the red channel of the OSC image for g_combined
            if self.ha_image is None:
                self.ha_image = r_channel
            # If OIII is None, use the green channel of the OSC image for b_combined
            if self.oiii_image is None:
                self.oiii_image = g_channel
            # If SII is None, use the red channel of the OSC image for r_combined
            if self.sii_image is None:
                self.sii_image = r_channel

            # Combined RGB channels with defaults as fallbacks
            r_combined = 0.5 * r_channel + 0.5 * self.sii_image
            g_combined = self.ha_to_oii_ratio * self.ha_image + (1 - self.ha_to_oii_ratio) * g_channel
            b_combined = b_channel
        else:
            # If no OSC image, use Ha, OIII, and SII images directly
            r_combined = 0.5 * self.ha_image + 0.5 * (self.sii_image if self.sii_image is not None else self.ha_image)
            g_combined = self.ha_to_oii_ratio * self.ha_image + (1 - self.ha_to_oii_ratio) * self.oiii_image
            b_combined = self.oiii_image

        # Stack the channels to create an RGB image
        combined_image = np.stack((r_combined, g_combined, b_combined), axis=-1)

        # Apply star stretch if enabled
        if self.enable_star_stretch:
            combined_image = self.apply_star_stretch(combined_image)

        # Apply SCNR (remove green cast)
        combined_image = self.apply_scnr(combined_image)
        self.preview_generated.emit(combined_image)


    def apply_star_stretch(self, image):
        stretched = ((3 ** self.stretch_factor) * image) / ((3 ** self.stretch_factor - 1) * image + 1)
        return np.clip(stretched, 0, 1)

    def apply_scnr(self, image):
        green_channel = image[..., 1]
        max_rg = np.maximum(image[..., 0], image[..., 2])
        green_channel[green_channel > max_rg] = max_rg[green_channel > max_rg]
        image[..., 1] = green_channel
        return image

class HaloBGonTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.image = None  # Selected image
        self.filename = None  # Store the selected file path
        self.preview_image = None  # Store the preview result
        self.zoom_factor = 0.25  # Initialize zoom factor for preview scaling
        self.dragging = False
        self.is_mono = True
        self.last_pos = None
        self.processing_thread = None  # For background processing
        self.original_header = None

    def initUI(self):
        main_layout = QHBoxLayout()

        # Left column for controls
        left_widget = QWidget(self)
        left_layout = QVBoxLayout(left_widget)
        left_widget.setFixedWidth(400)  # Fixed width for left column

        # Instructions label
        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Select a stars-only image.
            2. Adjust the reduction amount as needed.
            3. Click Execute to apply the halo reduction.
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # File selection button
        self.fileButton = QPushButton("Load Image", self)
        self.fileButton.clicked.connect(self.selectImage)
        left_layout.addWidget(self.fileButton)

        self.fileLabel = QLabel('No file selected', self)
        left_layout.addWidget(self.fileLabel)

        # Reduction amount slider
        self.reductionLabel = QLabel("Reduction Amount:", self)
        self.reductionSlider = QSlider(Qt.Horizontal, self)
        self.reductionSlider.setMinimum(0)
        self.reductionSlider.setMaximum(3)
        self.reductionSlider.setValue(1)
        self.reductionSlider.setToolTip("Adjust the amount of halo reduction (Extra Low, Low, Medium, High)")
        left_layout.addWidget(self.reductionLabel)
        left_layout.addWidget(self.reductionSlider)

        # Linear data checkbox
        self.linearDataCheckbox = QCheckBox("Linear Data", self)
        self.linearDataCheckbox.setToolTip("Check if the data is linear")
        left_layout.addWidget(self.linearDataCheckbox)

        # Progress indicator (spinner) label
        self.spinnerLabel = QLabel(self)
        self.spinnerLabel.setAlignment(Qt.AlignCenter)
        # Use the resource path function to access the GIF
        self.spinnerMovie = QMovie(resource_path("spinner.gif"))  # Updated path
        self.spinnerLabel.setMovie(self.spinnerMovie)
        self.spinnerLabel.hide()  # Hide spinner by default
        left_layout.addWidget(self.spinnerLabel)

        # Execute and Save buttons
        self.executeButton = QPushButton("Refresh Preview", self)
        self.executeButton.clicked.connect(self.generatePreview)
        left_layout.addWidget(self.executeButton)

        self.saveButton = QPushButton("Save Image", self)
        self.saveButton.clicked.connect(self.saveImage)
        left_layout.addWidget(self.saveButton)

        # Zoom in and out buttons
        zoom_layout = QHBoxLayout()
        self.zoomInButton = QPushButton("Zoom In", self)
        self.zoomInButton.clicked.connect(self.zoomIn)
        zoom_layout.addWidget(self.zoomInButton)

        self.zoomOutButton = QPushButton("Zoom Out", self)
        self.zoomOutButton.clicked.connect(self.zoomOut)
        zoom_layout.addWidget(self.zoomOutButton)
        left_layout.addLayout(zoom_layout)

                # Footer
        footer_label = QLabel("""
            Written by Franklin Marek<br>
            <a href='http://www.setiastro.com'>www.setiastro.com</a>
        """)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setOpenExternalLinks(True)
        footer_label.setStyleSheet("font-size: 10px;")
        left_layout.addWidget(footer_label)

        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))
        main_layout.addWidget(left_widget)

        

        # Right side for the preview inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.viewport().installEventFilter(self)

        # QLabel for the image preview
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400, 400)

        main_layout.addWidget(self.scrollArea)
        self.setLayout(main_layout)

    def zoomIn(self):
        self.zoom_factor *= 1.2  # Increase zoom by 20%
        self.updatePreview(self.image)

    def zoomOut(self):
        self.zoom_factor /= 1.2  # Decrease zoom by 20%
        self.updatePreview(self.image)
    

    def selectImage(self):
        selected_file, _ = QFileDialog.getOpenFileName(self, "Select Stars Only Image", "", "Images (*.png *.tif *.fits *.fit)")
        if selected_file:
            try:
                # Match the StarStretchTab loading method here
                self.image, self.original_header, _, self.is_mono = load_image(selected_file)  # Load image with header
                self.filename = selected_file 
                self.fileLabel.setText(os.path.basename(selected_file))
                self.generatePreview()  # Generate preview after loading
            except Exception as e:
                self.fileLabel.setText(f"Error: {str(e)}")
                print(f"Failed to load image: {e}")


    def applyHaloReduction(self):
        if self.image is None:
            print("No image selected.")
            return

        reduction_amount = self.reductionSlider.value()
        is_linear = self.linearDataCheckbox.isChecked()

        # Show spinner and start background processing
        self.showSpinner()
        self.processing_thread = QThread()
        self.processing_worker = self.HaloProcessingWorker(self.image, reduction_amount, is_linear)
        self.processing_worker.moveToThread(self.processing_thread)
        self.processing_worker.processing_complete.connect(self.updateImage)
        self.processing_thread.started.connect(self.processing_worker.process)
        self.processing_thread.start()

    def updatePreview(self, processed_image):
        # Update the preview once the processing thread emits the result
        preview_image = (processed_image * 255).astype(np.uint8)
        h, w = preview_image.shape[:2]
        if preview_image.ndim == 3:
            q_image = QImage(preview_image.data, w, h, 3 * w, QImage.Format_RGB888)
        else:
            q_image = QImage(preview_image.data, w, h, w, QImage.Format_Grayscale8)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())

        # Hide the spinner after processing is done
        self.hideSpinner()

    def saveImage(self):
        if self.image is not None:
            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename)
            default_save_name = os.path.splitext(base_name)[0] + '_reduced.tif'
            original_dir = os.path.dirname(self.filename)

            # Open the save file dialog
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name), 
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Call save_image with the necessary parameters
                        save_image(self.image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # For non-TIFF/FITS formats, save directly without bit depth selection
                    save_image(self.image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')



    def showSpinner(self):
        self.spinnerLabel.show()
        self.spinnerMovie.start()

    def hideSpinner(self):
        self.spinnerLabel.hide()
        self.spinnerMovie.stop()

    # Updated generatePreview method in HaloBGonTab to use HaloProcessingThread
    def generatePreview(self):
        if self.image is not None and self.image.size > 0:
            # Show spinner before starting processing
            self.showSpinner()

            # Start background processing with HaloProcessingThread
            self.processing_thread = HaloProcessingThread(self.image, self.reductionSlider.value(), self.linearDataCheckbox.isChecked())
            self.processing_thread.preview_generated.connect(self.updatePreview)
            self.processing_thread.start()

    def updatePreview(self, processed_image):
        # Update the preview once the processing thread emits the result
        preview_image = (processed_image * 255).astype(np.uint8)
        h, w = preview_image.shape[:2]
        if preview_image.ndim == 3:
            q_image = QImage(preview_image.data, w, h, 3 * w, QImage.Format_RGB888)
        else:
            q_image = QImage(preview_image.data, w, h, w, QImage.Format_Grayscale8)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())

        # Hide the spinner after processing is done
        self.hideSpinner()

    def eventFilter(self, source, event):
        if event.type() == event.MouseButtonPress and event.button() == Qt.LeftButton:
            self.dragging = True
            self.last_pos = event.pos()
        elif event.type() == event.MouseButtonRelease and event.button() == Qt.LeftButton:
            self.dragging = False
        elif event.type() == event.MouseMove and self.dragging:
            delta = event.pos() - self.last_pos
            self.scrollArea.horizontalScrollBar().setValue(self.scrollArea.horizontalScrollBar().value() - delta.x())
            self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().value() - delta.y())
            self.last_pos = event.pos()

        return super().eventFilter(source, event)


    def create_lightness_mask(image):
        # Convert to grayscale to get the lightness mask
        lightness_mask = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)
        
        # Apply Unsharp Mask
        blurred = cv2.GaussianBlur(lightness_mask, (0, 0), sigmaX=2)
        lightness_mask = cv2.addWeighted(lightness_mask, 1.66, blurred, -0.66, 0)
        
        # Normalize to the [0, 1] range
        return lightness_mask / 255.0

    def createDuplicateImage(self, original):
        return np.copy(original)

    def invert_mask(mask):
        return 1.0 - mask  # Assuming mask is normalized between 0 and 1


    def apply_mask_to_image(image, mask):
        # Ensure mask is 3-channel to match the image dimensions
        mask_rgb = np.stack([mask] * 3, axis=-1)
        return cv2.multiply(image, mask_rgb)


    def apply_curves_to_image(image, reduction_amount):
        # Define the curve based on reduction amount
        if reduction_amount == 0:
            curve = [int((i / 255.0) ** 0.575 * 255) for i in range(256)]
        else:
            curve = [int((i / 255.0) ** 0.4 * 255) for i in range(256)]
        
        lut = np.array(curve, dtype=np.uint8)
        return cv2.LUT((image * 255).astype(np.uint8), lut).astype(np.float32) / 255.0


    def load_image(self, filename):
        original_header = None
        file_extension = filename.split('.')[-1].lower()

        # Handle different file types and normalize them to [0, 1] range
        if file_extension in ['tif', 'tiff']:
            image = tiff.imread(filename).astype(np.float32) / 65535.0  # For 16-bit TIFF images
        elif file_extension == 'png':
            image = np.array(Image.open(filename).convert('RGB')).astype(np.float32) / 255.0  # Normalize to [0, 1]
        elif file_extension in ['fits', 'fit']:
            with fits.open(filename) as hdul:
                image = hdul[0].data.astype(np.float32)
                original_header = hdul[0].header
                # Normalize if data is 16-bit or higher
                if image.max() > 1:
                    image /= np.max(image)
        else:
            raise ValueError(f"Unsupported file format: {file_extension}")

        return image, original_header

    def save_image(self, image, filename, file_format, bit_depth="16-bit", original_header=None):
        img = Image.fromarray((image * 255).astype(np.uint8))
        img.save(filename)

class HaloProcessingThread(QThread):
    preview_generated = pyqtSignal(np.ndarray)

    def __init__(self, image, reduction_amount, is_linear):
        super().__init__()
        self.image = image
        self.reduction_amount = reduction_amount
        self.is_linear = is_linear

    def run(self):
        processed_image = self.applyHaloReduction(self.image, self.reduction_amount, self.is_linear)
        self.preview_generated.emit(processed_image)

    def applyHaloReduction(self, image, reduction_amount, is_linear):
        if is_linear:
            image = image ** (1 / 2.2)  # Gamma correction for linear data

        # Generate the lightness mask
        lightness_mask = self.createLightnessMask(image)
        inverted_mask = 1.0 - lightness_mask  # Invert the lightness mask

        # Duplicate and process the mask for halo reduction effect
        duplicated_mask = self.createDuplicateMask(lightness_mask)

        # Subtract duplicated mask from inverted mask to enhance halos
        enhanced_mask = inverted_mask - duplicated_mask * reduction_amount * 0.33  # Adjust factor based on reduction_amount

        # Apply the mask and curves transformation to the image
        masked_image = self.applyMaskToImage(image, enhanced_mask)
        final_image = self.applyCurvesToImage(masked_image, reduction_amount)

        if is_linear:
            final_image = final_image ** 2.2  # Restore gamma

        return np.clip(final_image, 0, 1)  # Ensure image stays within valid range

    def createLightnessMask(self, image):
        # Convert image to grayscale to create a lightness mask
        lightness_mask = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY).astype(np.float32) / 255.0
        blurred = cv2.GaussianBlur(lightness_mask, (0, 0), sigmaX=2)
        return cv2.addWeighted(lightness_mask, 1.66, blurred, -0.66, 0)

    def createDuplicateMask(self, mask):
        # Duplicate the mask and apply additional processing (simulating MMT)
        duplicated_mask = cv2.GaussianBlur(mask, (0, 0), sigmaX=2)
        return duplicated_mask

    def applyMaskToImage(self, image, mask):
        # Blend the original image with the mask based on the reduction level
        mask_rgb = np.stack([mask] * 3, axis=-1)  # Convert to 3-channel
        return cv2.multiply(image, mask_rgb)

    def applyCurvesToImage(self, image, reduction_amount):
        # Apply a curves transformation based on reduction_amount
        if reduction_amount == 0:
            # Extra Low setting, mild curve
            curve = [int((i / 255.0) ** 1.2 * 255) for i in range(256)]
        elif reduction_amount == 1:
            # Low setting, slightly stronger darkening
            curve = [int((i / 255.0) ** 1.5 * 255) for i in range(256)]
        elif reduction_amount == 2:
            # Medium setting, moderate darkening
            curve = [int((i / 255.0) ** 1.8 * 255) for i in range(256)]
        else:
            # High setting, strong darkening effect
            curve = [int((i / 255.0) ** 2.2 * 255) for i in range(256)]

        # Apply the curve transformation as a lookup table
        lut = np.array(curve, dtype=np.uint8)
        transformed_image = cv2.LUT((image * 255).astype(np.uint8), lut).astype(np.float32) / 255.0
        return transformed_image



class ContinuumSubtractTab(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.nb_image = None  # Changed from ha_image to nb_image
        self.filename = None  # Store the selected file path
        self.is_mono = True
        self.continuum_image = None  # Changed from red_continuum_image to continuum_image
        self.processing_thread = None  # For background processing
        self.combined_image = None  # Store the result of the continuum subtraction
        self.zoom_factor = 0.25  # Initial zoom factor

    def initUI(self):
        main_layout = QHBoxLayout()

        # Left side controls
        left_widget = QWidget(self)
        left_layout = QVBoxLayout(left_widget)
        left_widget.setFixedWidth(400)  # Fixed width for left column

        # Instruction box
        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Load your NB and Continuum images.
            2. Select for optional linear only output.
            3. Click Execute to perform continuum subtraction.
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # File Selection Buttons
        self.nb_button = QPushButton("Load NB Image")
        self.nb_button.clicked.connect(lambda: self.selectImage("nb"))
        self.nb_label = QLabel("No NB image selected")  # Updated label
        left_layout.addWidget(self.nb_button)
        left_layout.addWidget(self.nb_label)

        self.continuum_button = QPushButton("Load Continuum Image")
        self.continuum_button.clicked.connect(lambda: self.selectImage("continuum"))
        self.continuum_label = QLabel("No Continuum image selected")  # Updated label
        left_layout.addWidget(self.continuum_button)
        left_layout.addWidget(self.continuum_label)

        self.linear_output_checkbox = QCheckBox("Output Linear Image Only")
        left_layout.addWidget(self.linear_output_checkbox)

        # Progress indicator (spinner) label
        self.spinnerLabel = QLabel(self)
        self.spinnerLabel.setAlignment(Qt.AlignCenter)
        # Use the resource path function to access the GIF
        self.spinnerMovie = QMovie(resource_path("spinner.gif"))  # Updated path
        self.spinnerLabel.setMovie(self.spinnerMovie)
        self.spinnerLabel.hide()  # Hide spinner by default
        left_layout.addWidget(self.spinnerLabel)

        # Status label to show what is happening in the background
        self.statusLabel = QLabel(self)
        self.statusLabel.setAlignment(Qt.AlignCenter)
        left_layout.addWidget(self.statusLabel)


        # Execute Button
        self.execute_button = QPushButton("Execute")
        self.execute_button.clicked.connect(self.startContinuumSubtraction)
        left_layout.addWidget(self.execute_button)

        # Zoom In and Zoom Out Buttons
        zoom_layout = QHBoxLayout()
        self.zoomInButton = QPushButton("Zoom In")
        self.zoomInButton.clicked.connect(self.zoom_in)
        zoom_layout.addWidget(self.zoomInButton)

        self.zoomOutButton = QPushButton("Zoom Out")
        self.zoomOutButton.clicked.connect(self.zoom_out)
        zoom_layout.addWidget(self.zoomOutButton)
        left_layout.addLayout(zoom_layout)

        # Save Button
        self.save_button = QPushButton("Save Continuum Subtracted Image")
        self.save_button.clicked.connect(self.save_continuum_subtracted)
        left_layout.addWidget(self.save_button)

        # Footer
        footer_label = QLabel("""
            Written by Franklin Marek<br>
            <a href='http://www.setiastro.com'>www.setiastro.com</a>
        """)
        footer_label.setAlignment(Qt.AlignCenter)
        footer_label.setOpenExternalLinks(True)
        footer_label.setStyleSheet("font-size: 10px;")
        left_layout.addWidget(footer_label)

        # Spacer to push elements to the top
        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Add left widget layout to the main layout
        main_layout.addWidget(left_widget)

        # Right side for the preview inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.viewport().installEventFilter(self)

        # QLabel for the image preview
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400, 400)

        main_layout.addWidget(self.scrollArea)
        self.setLayout(main_layout)

    def selectImage(self, image_type):
        selected_file, _ = QFileDialog.getOpenFileName(self, "Select Image", "", "Images (*.png *.tif *.fits *.fit)")
        if selected_file:
            try:
                image, original_header, _, _ = load_image(selected_file)  # Load image with header
                self.filename = selected_file
                if image_type == "nb":
                    self.nb_image = image
                    self.nb_label.setText(os.path.basename(selected_file))  # Updated label
                elif image_type == "continuum":
                    self.continuum_image = image
                    self.continuum_label.setText(os.path.basename(selected_file))  # Updated label
            except Exception as e:
                print(f"Failed to load {image_type} image: {e}")
                if image_type == "nb":
                    self.nb_label.setText("Error loading NB image")
                elif image_type == "continuum":
                    self.continuum_label.setText("Error loading Continuum image")

    def startContinuumSubtraction(self):
        if self.nb_image is not None and self.continuum_image is not None:
            # Show spinner and start background processing
            self.showSpinner()
            self.processing_thread = ContinuumProcessingThread(self.nb_image, self.continuum_image,
                                                            self.linear_output_checkbox.isChecked())
            self.processing_thread.processing_complete.connect(self.display_image)
            self.processing_thread.finished.connect(self.hideSpinner)
            self.processing_thread.status_update.connect(self.update_status_label)
            self.processing_thread.start()
        else:
            print("Please select both NB and Continuum images.")

    def update_status_label(self, message):
        self.statusLabel.setText(message)

    def zoom_in(self):
        self.zoom_factor *= 1.2
        self.update_preview()

    def zoom_out(self):
        self.zoom_factor /= 1.2
        self.update_preview()

    def update_preview(self):
        if self.combined_image is not None:
            self.display_image(self.combined_image)        

    def load_image(self, filename):
        # Placeholder for actual image loading logic
        image = cv2.imread(filename, cv2.IMREAD_GRAYSCALE).astype(np.float32) / 255.0
        return image, None, None, None
    
    def save_continuum_subtracted(self):
        if self.combined_image is not None:
            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename)
            default_save_name = os.path.splitext(base_name)[0] + '_continuumsubtracted.tif'
            original_dir = os.path.dirname(self.filename)

            # Open the save file dialog
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name), 
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Call save_image with the necessary parameters
                        save_image(self.image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # For non-TIFF/FITS formats, save directly without bit depth selection
                    save_image(self.image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')



    def display_image(self, processed_image):
        self.combined_image = processed_image

        # Convert the processed image to a displayable format
        preview_image = (processed_image * 255).astype(np.uint8)
        
        # Check if the image is mono or RGB
        if preview_image.ndim == 2:  # Mono image
            # Create a 3-channel RGB image by duplicating the single channel
            preview_image = np.stack([preview_image] * 3, axis=-1)  # Stack to create RGB

        h, w = preview_image.shape[:2]

        # Change the format to RGB888 for displaying an RGB image
        q_image = QImage(preview_image.data, w, h, 3 * w, QImage.Format_RGB888)

        pixmap = QPixmap.fromImage(q_image)
        scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        self.imageLabel.setPixmap(scaled_pixmap)
        self.imageLabel.resize(scaled_pixmap.size())


    def showSpinner(self):
        self.spinnerLabel.show()
        self.spinnerMovie.start()

    def hideSpinner(self):
        self.spinnerLabel.hide()
        self.spinnerMovie.stop()

    def eventFilter(self, source, event):
        if source is self.scrollArea.viewport():
            if event.type() == event.MouseButtonPress and event.button() == Qt.LeftButton:
                self.dragging = True
                self.last_pos = event.pos()
            elif event.type() == event.MouseButtonRelease and event.button() == Qt.LeftButton:
                self.dragging = False
            elif event.type() == event.MouseMove and self.dragging:
                delta = event.pos() - self.last_pos
                self.scrollArea.horizontalScrollBar().setValue(self.scrollArea.horizontalScrollBar().value() - delta.x())
                self.scrollArea.verticalScrollBar().setValue(self.scrollArea.verticalScrollBar().value() - delta.y())
                self.last_pos = event.pos()

        return super().eventFilter(source, event)


class ContinuumProcessingThread(QThread):
    processing_complete = pyqtSignal(np.ndarray)
    status_update = pyqtSignal(str)

    def __init__(self, nb_image, continuum_image, output_linear):
        super().__init__()
        self.nb_image = nb_image
        self.continuum_image = continuum_image
        self.output_linear = output_linear
        self.background_reference = None  # Store the background reference



    def run(self):
        # Ensure both images are mono
        if self.nb_image.ndim == 3 and self.nb_image.shape[2] == 3:
            self.nb_image = self.nb_image[..., 0]  # Take one channel for the NB image

        if self.continuum_image.ndim == 3 and self.continuum_image.shape[2] == 3:
            self.continuum_image = self.continuum_image[..., 0]  # Take one channel for the continuum image

        # Create RGB image
        r_combined = self.nb_image  # Use the normalized NB image as the Red channel
        g_combined = self.continuum_image # Use the normalized continuum image as the Green channel
        b_combined = self.continuum_image  # Use the normalized continuum image as the Blue channel


        # Stack the channels into a single RGB image
        combined_image = np.stack((r_combined, g_combined, b_combined), axis=-1)

        self.status_update.emit("Performing background neutralization...")
        QCoreApplication.processEvents()
            # Perform background neutralization
        self.background_neutralization(combined_image)

        # Normalize the red channel to the green channel
        combined_image[..., 0] = self.normalize_channel(combined_image[..., 0], combined_image[..., 1])

        # Perform continuum subtraction
        linear_image = combined_image[..., 0] - 0.9*(combined_image[..., 1]-np.median(combined_image[..., 1]))

            # Check if the Output Linear checkbox is checked
        if self.output_linear:
            # Emit the linear image for preview
            self.processing_complete.emit(np.clip(linear_image, 0, 1))
            return  # Exit the method if we only want to output the linear image

        self.status_update.emit("Subtraction complete.")
        QCoreApplication.processEvents()

        # Perform statistical stretch
        target_median = 0.25
        stretched_image = stretch_color_image(linear_image, target_median, True, False)

        # Final image adjustment
        final_image = stretched_image - 0.7*np.median(stretched_image)

        # Clip the final image to stay within [0, 1]
        final_image = np.clip(final_image, 0, 1)

        # Applies Curves Boost
        final_image = apply_curves_adjustment(final_image, np.median(final_image), 0.5)

        self.status_update.emit("Linear to Non-Linear Stretch complete.")
        QCoreApplication.processEvents()
        # Emit the final image for preview
        self.processing_complete.emit(final_image)

    def background_neutralization(self, rgb_image):
        height, width, _ = rgb_image.shape
        num_boxes = 200
        box_size = 25
        iterations = 25

        boxes = [(np.random.randint(0, height - box_size), np.random.randint(0, width - box_size)) for _ in range(num_boxes)]
        best_means = np.full(num_boxes, np.inf)

        for _ in range(iterations):
            for i, (y, x) in enumerate(boxes):
                if y + box_size <= height and x + box_size <= width:
                    patch = rgb_image[y:y + box_size, x:x + box_size]
                    patch_median = np.median(patch) if patch.size > 0 else np.inf

                    if patch_median < best_means[i]:
                        best_means[i] = patch_median

                    surrounding_values = []
                    for dy in [-1, 0, 1]:
                        for dx in [-1, 0, 1]:
                            surrounding_y = y + dy * box_size
                            surrounding_x = x + dx * box_size
                            
                            if (0 <= surrounding_y < height - box_size) and (0 <= surrounding_x < width - box_size):
                                surrounding_patch = rgb_image[surrounding_y:surrounding_y + box_size, surrounding_x:surrounding_x + box_size]
                                if surrounding_patch.size > 0:
                                    surrounding_values.append(np.median(surrounding_patch))

                    if surrounding_values:
                        dimmest_index = np.argmin(surrounding_values)
                        new_y = y + (dimmest_index // 3 - 1) * box_size
                        new_x = x + (dimmest_index % 3 - 1) * box_size
                        boxes[i] = (new_y, new_x)

        # After iterations, find the darkest box median
        darkest_value = np.inf
        background_box = None

        for box in boxes:
            y, x = box
            if y + box_size <= height and x + box_size <= width:
                patch = rgb_image[y:y + box_size, x:y + box_size]
                patch_median = np.median(patch) if patch.size > 0 else np.inf

                if patch_median < darkest_value:
                    darkest_value = patch_median
                    background_box = patch

        if background_box is not None:
            self.background_reference = np.median(background_box.reshape(-1, 3), axis=0)
            
            # Adjust the channels based on the median reference
            channel_medians = np.median(rgb_image, axis=(0, 1))

            # Adjust channels based on the red channel
            for channel in range(3):
                if self.background_reference[channel] < channel_medians[channel]:
                    pedestal = channel_medians[channel] - self.background_reference[channel]
                    rgb_image[..., channel] += pedestal

            # Specifically adjust G and B to match R
            r_median = self.background_reference[0]
            for channel in [1, 2]:  # Green and Blue channels
                if self.background_reference[channel] < r_median:
                    rgb_image[..., channel] += (r_median - self.background_reference[channel])

        self.status_update.emit("Background neutralization complete.")
        QCoreApplication.processEvents()
        return rgb_image
    
    def normalize_channel(self, image_channel, reference_channel):
        mad_image = np.mean(np.abs(image_channel - np.mean(image_channel)))
        mad_ref = np.mean(np.abs(reference_channel - np.mean(reference_channel)))

        median_image = np.median(image_channel)
        median_ref = np.median(reference_channel)

        # Apply the normalization formula
        normalized_channel = (
            image_channel * mad_ref / mad_image
            - (mad_ref / mad_image) * median_image
            + median_ref
        )

        self.status_update.emit("Color calibration complete.")
        QCoreApplication.processEvents()
        return np.clip(normalized_channel, 0, 1)  



    def continuum_subtraction(self, rgb_image):
        red_channel = rgb_image[..., 0]
        green_channel = rgb_image[..., 1]
        
        # Determine Q based on the selection (modify condition based on actual UI element)
        Q = 0.9 if self.output_linear else 1.0

        # Perform the continuum subtraction
        median_green = np.median(green_channel)
        result_image = red_channel - Q * (green_channel - median_green)
        
        return np.clip(result_image, 0, 1)  # Ensure values stay within [0, 1]





# Function to load and convert images to 32-bit floating point
def load_image(filename):
    bit_depth = None  # Initialize bit depth to None
    is_mono = True  # Assume monochrome by default    
    original_header = None  # Initialize an empty header for FITS files
    if filename.lower().endswith('.png'):
        img = Image.open(filename).convert('RGB')  # Ensures it's RGB
        img_array = np.array(img, dtype=np.float32) / 255.0  # Normalize to [0, 1]
    elif filename.lower().endswith('.tiff') or filename.lower().endswith('.tif'):
        img_array = tiff.imread(filename)
        # Handle different bit depths
        if img_array.dtype == np.uint8:
            img_array = img_array.astype(np.float32) / 255.0  # Normalize 8-bit to [0, 1]
        elif img_array.dtype == np.uint16:
            img_array = img_array.astype(np.float32) / 65535.0  # Normalize 16-bit to [0, 1]
        elif img_array.dtype == np.uint32:
            img_array = img_array.astype(np.float32) / 4294967295.0  # Normalize 32-bit unsigned integer to [0, 1]
        elif img_array.dtype == np.float32:
            img_array = img_array  # Already in 32-bit floating point, no need to convert
        else:
            raise ValueError("Unsupported TIFF format!")
    elif filename.lower().endswith('.fits') or filename.lower().endswith('.fit'):
        with fits.open(filename) as hdul:
            img_array = hdul[0].data
            original_header = hdul[0].header  # Capture the FITS header
            bit_depth = None

            # Determine bit depth and apply necessary transformations
            if img_array.dtype == np.uint16:
                bit_depth = "16-bit"
                img_array = img_array.astype(np.float32) / 65535.0  # Normalize 16-bit to [0, 1]
            elif img_array.dtype == np.uint32:
                bit_depth = "32-bit unsigned"
                bzero = original_header.get('BZERO', 0)
                bscale = original_header.get('BSCALE', 1)
                img_array = img_array.astype(np.float32) * bscale + bzero

                # Normalize to [0, 1] based on range
                image_min = img_array.min()
                image_max = img_array.max()
                img_array = (img_array - image_min) / (image_max - image_min)
            elif img_array.dtype == np.float32:
                bit_depth = "32-bit floating point"
                # No normalization needed for 32-bit float

            # Handle 3D FITS data (e.g., RGB or multi-layered data)
            if img_array.ndim == 3 and img_array.shape[0] == 3:
                img_array = np.transpose(img_array, (1, 2, 0))  # Reorder to (height, width, channels)
                is_mono = False
            elif img_array.ndim == 2:
                img_array = np.stack([img_array] * 3, axis=-1)  # Convert grayscale to 3-channel for consistency
                is_mono = True
            else:
                raise ValueError("Unsupported FITS format!")
    else:
        raise ValueError("Unsupported file format!")

    return img_array, original_header, bit_depth, is_mono  # Return the image array, header, bit depth, and is_mono flag






def save_image(img_array, filename, original_format, bit_depth=None, original_header=None, is_mono=False):
    img_array = ensure_native_byte_order(img_array)  # Apply native byte order correction if needed

    if original_format == 'png':
        img = Image.fromarray((img_array * 255).astype(np.uint8))  # Convert to 8-bit and save as PNG
        img.save(filename)
    elif original_format in ['tiff', 'tif']:
        if bit_depth == "16-bit":
            tiff.imwrite(filename, (img_array * 65535).astype(np.uint16))  # Save as 16-bit TIFF
        elif bit_depth == "32-bit unsigned":
            tiff.imwrite(filename, (img_array * 4294967295).astype(np.uint32))  # Save as 32-bit unsigned TIFF
        elif bit_depth == "32-bit floating point":
            tiff.imwrite(filename, img_array.astype(np.float32))  # Save as 32-bit floating point TIFF
    elif original_format in ['fits', 'fit']:
        # For grayscale (mono) FITS images
        if is_mono:
            if bit_depth == "16-bit":
                img_array_fits = (img_array[:, :, 0] * 65535).astype(np.uint16)
            elif bit_depth == "32-bit unsigned":
                img_array_fits = (img_array[:, :, 0] * 4294967295).astype(np.uint32)
            elif bit_depth == "32-bit floating point":
                img_array_fits = img_array[:, :, 0].astype(np.float32)
            hdu = fits.PrimaryHDU(img_array_fits, header=original_header)
        else:
            # Transpose RGB image to (channels, height, width) for FITS format
            img_array_fits = np.transpose(img_array, (2, 0, 1))
            if bit_depth == "16-bit":
                img_array_fits = (img_array_fits * 65535).astype(np.uint16)
            elif bit_depth == "32-bit unsigned":
                img_array_fits = (img_array_fits * 4294967295).astype(np.uint32)
            elif bit_depth == "32-bit floating point":
                img_array_fits = img_array_fits.astype(np.float32)

            # Update the original header with correct dimensions for multi-channel images
            original_header['NAXIS'] = 3
            original_header['NAXIS1'] = img_array_fits.shape[2]  # Width
            original_header['NAXIS2'] = img_array_fits.shape[1]  # Height
            original_header['NAXIS3'] = img_array_fits.shape[0]  # Channels

            hdu = fits.PrimaryHDU(img_array_fits, header=original_header)

        hdu.writeto(filename, overwrite=True)
    else:
        raise ValueError("Unsupported file format!")




def stretch_mono_image(image, target_median, normalize=False, apply_curves=False, curves_boost=0.0):
    black_point = max(np.min(image), np.median(image) - 2.7 * np.std(image))
    rescaled_image = (image - black_point) / (1 - black_point)
    median_image = np.median(rescaled_image)
    stretched_image = ((median_image - 1) * target_median * rescaled_image) / (median_image * (target_median + rescaled_image - 1) - target_median * rescaled_image)
    
    if apply_curves:
        stretched_image = apply_curves_adjustment(stretched_image, target_median, curves_boost)
    
    if normalize:
        stretched_image = stretched_image / np.max(stretched_image)
    
    return np.clip(stretched_image, 0, 1)


def stretch_color_image(image, target_median, linked=True, normalize=False, apply_curves=False, curves_boost=0.0):
    if linked:
        combined_median = np.median(image)
        combined_std = np.std(image)
        black_point = max(np.min(image), combined_median - 2.7 * combined_std)
        rescaled_image = (image - black_point) / (1 - black_point)
        median_image = np.median(rescaled_image)
        stretched_image = ((median_image - 1) * target_median * rescaled_image) / (median_image * (target_median + rescaled_image - 1) - target_median * rescaled_image)
    else:
        stretched_image = np.zeros_like(image)
        for channel in range(image.shape[2]):
            black_point = max(np.min(image[..., channel]), np.median(image[..., channel]) - 2.7 * np.std(image[..., channel]))
            rescaled_channel = (image[..., channel] - black_point) / (1 - black_point)
            median_channel = np.median(rescaled_channel)
            stretched_image[..., channel] = ((median_channel - 1) * target_median * rescaled_channel) / (median_channel * (target_median + rescaled_channel - 1) - target_median * rescaled_channel)
    
    if apply_curves:
        stretched_image = apply_curves_adjustment(stretched_image, target_median, curves_boost)
    
    if normalize:
        stretched_image = stretched_image / np.max(stretched_image)
    
    return np.clip(stretched_image, 0, 1)


def apply_curves_adjustment(image, target_median, curves_boost):
    curve = [
        [0.0, 0.0],
        [0.5 * target_median, 0.5 * target_median],
        [target_median, target_median],
        [(1 / 4 * (1 - target_median) + target_median), 
         np.power((1 / 4 * (1 - target_median) + target_median), (1 - curves_boost))],
        [(3 / 4 * (1 - target_median) + target_median), 
         np.power(np.power((3 / 4 * (1 - target_median) + target_median), (1 - curves_boost)), (1 - curves_boost))],
        [1.0, 1.0]
    ]
    adjusted_image = np.interp(image, [p[0] for p in curve], [p[1] for p in curve])
    return adjusted_image

def resource_path(relative_path):
    """ Get the absolute path to the resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temporary folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def ensure_native_byte_order(array):
    """
    Ensures that the array is in the native byte order.
    If the array is in a non-native byte order, it will convert it.
    """
    if array.dtype.byteorder == '=':  # Already in native byte order
        return array
    elif array.dtype.byteorder in ('<', '>'):  # Non-native byte order
        return array.byteswap().newbyteorder()
    return array


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = AstroEditingSuite()
    window.show()
    sys.exit(app.exec_())
