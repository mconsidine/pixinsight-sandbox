import numpy as np
import os
import tifffile as tiff
from PIL import Image
from astropy.io import fits
import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QPushButton, QLabel, QSlider, QCheckBox, 
                             QFileDialog, QVBoxLayout, QHBoxLayout, QSpacerItem, QSizePolicy)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QImage
from PyQt5.QtWidgets import QScrollArea
from PyQt5.QtCore import Qt, QPoint
from PyQt5.QtWidgets import QInputDialog


# SetiAstro Copyright Notice
def display_setiastro_copyright():
    print((r"""
 *#        ___     __      ___       __                              #
 *#       / __/___/ /__   / _ | ___ / /________                      #
 *#      _\ \/ -_) _ _   / __ |(_-</ __/ __/ _ \                     #
 *#     /___/\__/_//_/  /_/ |_/___/\__/__/ \___/                     #
 *#                                                                  #
 *#              Statistical Stretch - V1.4                          #
 *#                                                                  #
 *#                         SetiAstro                                #
 *#                    Copyright © 2024                              #
 *#                                                                  #
    """))

# Display the copyright notice at the start of the script
display_setiastro_copyright()

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
        # Save as FITS file with header information
        if is_mono:
            # For grayscale, save only the first channel with header information
            if bit_depth == "16-bit":
                img_array_fits = (img_array[:, :, 0] * 65535).astype(np.uint16)
            elif bit_depth == "32-bit unsigned":
                img_array_fits = img_array[:, :, 0].astype(np.float32)
            elif bit_depth == "32-bit floating point":
                img_array_fits = img_array[:, :, 0].astype(np.float32)
            hdu = fits.PrimaryHDU(img_array_fits, header=original_header)
        else:
            # Transpose RGB image back to (channels, height, width) for FITS saving
            img_array_fits = np.transpose(img_array, (2, 0, 1))

            if bit_depth == "16-bit":
                img_array_fits = (img_array_fits * 65535).astype(np.uint16)
            elif bit_depth == "32-bit unsigned":
                img_array_fits = img_array_fits.astype(np.float32)
                original_header['BITPIX'] = -32
            elif bit_depth == "32-bit floating point":
                img_array_fits = img_array_fits.astype(np.float32)

            # Update the original header to reflect the correct dimensions
            original_header['NAXIS'] = 3  # Number of axes
            original_header['NAXIS1'] = img_array_fits.shape[2]  # Width
            original_header['NAXIS2'] = img_array_fits.shape[1]  # Height
            original_header['NAXIS3'] = img_array_fits.shape[0]  # Number of channels

            hdu = fits.PrimaryHDU(img_array_fits, header=original_header)

        hdu.writeto(filename, overwrite=True)
    else:
        raise ValueError("Unsupported file format!")


# Curves adjustment function
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


# Mono image stretch function
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


# Color image stretch function
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


# PyQt GUI for the Statistical Stretch application
class ImageStretchApp(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.image = None
        self.filename = None
        self.zoom_factor = 1.0

    def initUI(self):
        self.setWindowTitle('Statistical Stretch - V1.4')
        main_layout = QHBoxLayout()  # Main layout is horizontal to allow for left and right sections

        # Left column (fixed width layout)
        left_widget = QWidget(self)  # Wrap left layout in a QWidget
        left_layout = QVBoxLayout(left_widget)

        # Set fixed width for the left column
        left_widget.setFixedWidth(400)  # You can adjust this width as needed

        # Instruction box at the top of the window
        instruction_box = QLabel(self)
        instruction_box.setText("""
            Instructions:
            1. Select an image to stretch.
            2. Adjust the target median and optional settings.
                0.10 is a good start for distinct objects eg. galaxies and pn
                0.25 is a good start for distended objects like nebula filling the image                               
            3. Preview the result.
            4. Save the stretched image in your desired format.

            Options:
            - Linked Stretch: Apply stretch to all channels together.
            - Normalize Image: Adjust the image range to [0, 1].
            - Apply Curves Adjustment: Add a final curves boost.

            SetiAstro Copyright © 2024
        """)
        instruction_box.setWordWrap(True)
        left_layout.addWidget(instruction_box)

        # File selection button
        self.fileButton = QPushButton('Select Image', self)
        self.fileButton.clicked.connect(self.openFileDialog)
        left_layout.addWidget(self.fileButton)

        # Label to show selected file name
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
        self.curvesCheckBox.stateChanged.connect(self.toggleCurvesSlider)  # Connect to toggle function
        left_layout.addWidget(self.curvesCheckBox)

        # Curves Boost slider (initially hidden)
        self.curvesBoostLabel = QLabel('Curves Boost: 0.00', self)
        self.curvesBoostSlider = QSlider(Qt.Horizontal)
        self.curvesBoostSlider.setMinimum(0)
        self.curvesBoostSlider.setMaximum(50)
        self.curvesBoostSlider.setValue(0)
        self.curvesBoostSlider.valueChanged.connect(self.updateCurvesBoostLabel)

        # Initially hide the curves boost slider and label
        self.curvesBoostLabel.hide()
        self.curvesBoostSlider.hide()

        left_layout.addWidget(self.curvesBoostLabel)
        left_layout.addWidget(self.curvesBoostSlider)

        # Preview button
        self.previewButton = QPushButton('Preview Stretch', self)
        self.previewButton.clicked.connect(self.previewStretch)
        left_layout.addWidget(self.previewButton)

        # Zoom buttons for preview
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

        left_layout.addSpacerItem(QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding))  # Add space to fill remaining

        # Add the left_widget (with fixed width) to the main layout
        main_layout.addWidget(left_widget)

        # Right side for the preview (stretchable) inside a QScrollArea
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        self.scrollArea.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        # QLabel inside QScrollArea for the image preview
        self.imageLabel = QLabel(self)
        self.imageLabel.setAlignment(Qt.AlignCenter)  # Center the image inside the scroll area
        self.scrollArea.setWidget(self.imageLabel)
        self.scrollArea.setMinimumSize(400,400)

        # Add the scroll area to the main layout so it stretches
        main_layout.addWidget(self.scrollArea)

        # Set the main layout for the entire window
        self.setLayout(main_layout)

        # Start zoom at 25%
        self.zoom_factor = 0.25

        # Enable mouse dragging
        self.scrollArea.viewport().setMouseTracking(True)
        self.scrollArea.viewport().installEventFilter(self)

        self.dragging = False
        self.last_pos = QPoint()

    # Override eventFilter to handle click and drag
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


    # Function to zoom in
    def zoom_in(self):
        self.zoom_factor *= 1.2
        self.update_preview()

    # Function to zoom out
    def zoom_out(self):
        self.zoom_factor /= 1.2
        self.update_preview()

    # Function to toggle the visibility of curves boost slider
    def toggleCurvesSlider(self, state):
        if state == Qt.Checked:
            self.curvesBoostLabel.show()
            self.curvesBoostSlider.show()
        else:
            self.curvesBoostLabel.hide()
            self.curvesBoostSlider.hide()

    def openFileDialog(self):
        # Updated file filter to include .fit along with .fits, .png, and .tiff
        self.filename, _ = QFileDialog.getOpenFileName(self, 'Open Image File', '', 
                                                    'Images (*.png *.tiff *.tif *.fits *.fit);;All Files (*)')
        if self.filename:
            self.fileLabel.setText(self.filename)
            self.image, self.original_header, self.bit_depth, self.is_mono = load_image(self.filename)  # Unpack all four values



    def updateMedianLabel(self, value):
        self.medianLabel.setText(f'Target Median: {value / 100:.2f}')

    def updateCurvesBoostLabel(self, value):
        self.curvesBoostLabel.setText(f'Curves Boost: {value / 100:.2f}')

    def previewStretch(self):
        if self.image is not None:
            target_median = self.medianSlider.value() / 100.0
            linked = self.linkedCheckBox.isChecked()
            normalize = self.normalizeCheckBox.isChecked()
            apply_curves = self.curvesCheckBox.isChecked()  # Check if curves should be applied
            curves_boost = self.curvesBoostSlider.value() / 100.0  # Get the curves boost value

            # Apply the stretch
            if self.image.ndim == 2:  # Mono image
                stretched_image = stretch_mono_image(self.image, target_median, normalize, apply_curves, curves_boost)
            else:  # Color image
                stretched_image = stretch_color_image(self.image, target_median, linked, normalize, apply_curves, curves_boost)

            # Store the stretched image for zooming
            self.stretched_image = stretched_image
            
            # Only set zoom_factor to 0.25 on the first preview
            if not hasattr(self, 'zoom_initialized') or not self.zoom_initialized:
                self.zoom_factor = 0.25
                self.zoom_initialized = True  # Mark that zoom has been initialized

            # Update the preview with the current zoom level
            self.update_preview()



    def update_preview(self):
        if hasattr(self, 'stretched_image'):
            img = (self.stretched_image * 255).astype(np.uint8)  # Scale to 8-bit values
            h, w = img.shape[:2]

            # Handle mono vs color images
            if img.ndim == 3:  # RGB image
                bytes_per_line = 3 * w
                q_image = QImage(img.tobytes(), w, h, bytes_per_line, QImage.Format_RGB888)
            else:  # Grayscale image
                bytes_per_line = w
                q_image = QImage(img.tobytes(), w, h, bytes_per_line, QImage.Format_Grayscale8)

            # Apply zoom factor
            pixmap = QPixmap.fromImage(q_image)
            scaled_pixmap = pixmap.scaled(pixmap.size() * self.zoom_factor, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            
            # Set scaled image to QLabel and resize QLabel
            self.imageLabel.setPixmap(scaled_pixmap)
            self.imageLabel.resize(scaled_pixmap.size())

            # Ensure the scroll area shows the correct area
            self.scrollArea.setWidget(self.imageLabel)



    def saveImage(self):
        if self.image is not None:
            target_median = self.medianSlider.value() / 100.0
            linked = self.linkedCheckBox.isChecked()
            normalize = self.normalizeCheckBox.isChecked()
            apply_curves = self.curvesCheckBox.isChecked()
            curves_boost = self.curvesBoostSlider.value() / 100.0

            # Perform the stretch based on the image type
            if self.image.ndim == 2:  # Mono image
                stretched_image = stretch_mono_image(self.image, target_median, normalize, apply_curves, curves_boost)
            else:  # Color image
                stretched_image = stretch_color_image(self.image, target_median, linked, normalize, apply_curves, curves_boost)

            # Pre-populate the save dialog with the original image name
            base_name = os.path.basename(self.filename)  # Extract the original image name
            default_save_name = os.path.splitext(base_name)[0] + '_stretched.tif'  # Set the default extension to .tif

            # Get the directory of the original image
            original_dir = os.path.dirname(self.filename)

            # Open file save dialog with default name and directory
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 
                'Save Image As', 
                os.path.join(original_dir, default_save_name),  # Pre-populate the directory and name
                'Images (*.tiff *.tif *.png *.fit *.fits);;All Files (*)'
            )

            # If the user selected a file, proceed with saving
            if save_filename:
                original_format = save_filename.split('.')[-1].lower()

                # For TIFF and FITS files, prompt the user to select the bit depth
                if original_format in ['tiff', 'tif', 'fits', 'fit']:
                    bit_depth_options = ["16-bit", "32-bit unsigned", "32-bit floating point"]
                    bit_depth, ok = QInputDialog.getItem(self, "Select Bit Depth", "Choose bit depth for saving:", bit_depth_options, 0, False)
                    
                    if ok and bit_depth:
                        # Pass the selected bit depth to the save function
                        save_image(stretched_image, save_filename, original_format, bit_depth, self.original_header, self.is_mono)
                        self.fileLabel.setText(f'Image saved as: {save_filename}')
                    else:
                        self.fileLabel.setText('Save canceled.')
                else:
                    # If it's not TIFF or FITS, just save the image directly
                    save_image(stretched_image, save_filename, original_format)
                    self.fileLabel.setText(f'Image saved as: {save_filename}')
            else:
                self.fileLabel.setText('Save canceled.')


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = ImageStretchApp()
    ex.show()
    sys.exit(app.exec_())
