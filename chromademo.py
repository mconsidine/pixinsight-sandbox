import cv2
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from threading import Timer


# Global variable to track if a slider is being moved
slider_active = False
text_original, text_saturation, text_chroma = None, None, None

def on_slider_release(event):
    """Update images only when releasing the slider."""
    global slider_active
    if slider_active:
        update_images()
        slider_active = False

def on_slider_change(val):
    """Set slider_active flag to True when moving the slider."""
    global slider_active
    slider_active = True

def adjust_saturation_hsv(image, factor):
    """Adjust saturation using HSV color space."""
    hsv = cv2.cvtColor(image, cv2.COLOR_RGB2HSV)  # Convert to HSV
    hsv[:, :, 1] = np.clip(hsv[:, :, 1] * factor, 0, 255)  # Scale saturation channel
    return cv2.cvtColor(hsv, cv2.COLOR_HSV2RGB)  # Convert back to RGB

def adjust_chroma_lab(image, chroma_factor=1.0):
    """Convert to LAB, scale chroma (a and b channels), and convert back to RGB."""
    # Convert to LAB color space
    lab = cv2.cvtColor(image, cv2.COLOR_RGB2LAB).astype(np.float32)
    
    # Scale L channel properly
    lab[..., 0] = lab[..., 0] * (100 / 255.0)
    
    # Adjust chroma (a and b channels)
    lab[..., 1] = (lab[..., 1] - 128) * chroma_factor + 128  # Scale chroma
    lab[..., 2] = (lab[..., 2] - 128) * chroma_factor + 128
    
    # Scale L channel back
    lab[..., 0] = lab[..., 0] * (255.0 / 100.0)
    
    # Clip values and convert back to uint8
    lab = np.clip(lab, 0, 255).astype(np.uint8)
    
    # Convert back to RGB
    return cv2.cvtColor(lab, cv2.COLOR_LAB2RGB)

# Function to load an image based on user selection
def load_image(image_choice):
    global image
    image_path = resource_path(image_choice)
    image = cv2.imread(image_path)
    if image is None:
        raise FileNotFoundError(f"Image file '{image_path}' not found!")
    image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

def update_bar_graph():
    """Update the RGB bar graphs and display combined brightness values."""
    global sample_x, sample_y, text_original, text_saturation, text_chroma

    # Sample point RGB values
    sample_rgb = image[sample_y, sample_x].astype(np.float32) / 255.0
    sample_rgb_saturated = sat_image[sample_y, sample_x].astype(np.float32) / 255.0
    sample_rgb_chroma = chroma_image[sample_y, sample_x].astype(np.float32) / 255.0

    # Combined brightness (mean of R, G, and B)
    brightness_original = np.mean(sample_rgb)
    brightness_saturated = np.mean(sample_rgb_saturated)
    brightness_chroma = np.mean(sample_rgb_chroma)

    # Update the bar graphs
    for bars, values in zip(
        [bars_original, bars_saturation, bars_chroma],
        [sample_rgb, sample_rgb_saturated, sample_rgb_chroma]
    ):
        for bar, value in zip(bars, values):
            bar.set_height(value)

    # Clear previous text
    if text_original: text_original.remove()
    if text_saturation: text_saturation.remove()
    if text_chroma: text_chroma.remove()

    # Add new text objects
    text_original = bar_ax_original.text(0.5, -0.3, f"Avg Brightness: {brightness_original:.2f}",
                                         transform=bar_ax_original.transAxes, ha='center', fontsize=10)
    text_saturation = bar_ax_saturation.text(0.5, -0.3, f"Avg Brightness: {brightness_saturated:.2f}",
                                             transform=bar_ax_saturation.transAxes, ha='center', fontsize=10)
    text_chroma = bar_ax_chroma.text(0.5, -0.3, f"Avg Brightness: {brightness_chroma:.2f}",
                                     transform=bar_ax_chroma.transAxes, ha='center', fontsize=10)

    # Redraw the graphs
    for bar_ax in [bar_ax_original, bar_ax_saturation, bar_ax_chroma]:
        bar_ax.relim()
        bar_ax.autoscale_view()

    fig.canvas.draw_idle()

def onclick(event):
    """Update the sample point based on mouse click."""
    global sample_x, sample_y
    if event.inaxes == ax_original:
        x, y = int(event.xdata), int(event.ydata)
        # Ensure the click is within image bounds
        if 0 <= x < image.shape[1] and 0 <= y < image.shape[0]:
            sample_x, sample_y = x, y
            update_bar_graph()

# Debounce mechanism for sliders
def debounce_update_images():
    global debounce_timer
    if debounce_timer is not None:
        debounce_timer.cancel()
    debounce_timer = Timer(0.5, update_images)
    debounce_timer.start()

def update_images():
    global sat_image, chroma_image

    # Update saturation image
    sat_image = adjust_saturation_hsv(image, sat_slider.val)
    ax_saturation.imshow(sat_image)

    # Update chroma image
    chroma_image = adjust_chroma_lab(image, chroma_slider.val)
    ax_chroma.imshow(chroma_image)

    update_bar_graph()  # Update graphs after image update

    fig.canvas.draw_idle()

def resource_path(relative_path):
    """Get the absolute path to a resource, works for PyInstaller bundled apps."""
    if hasattr(sys, '_MEIPASS'):
        # PyInstaller creates a temporary folder for bundled files
        return os.path.join(sys._MEIPASS, relative_path)
    return os.path.join(os.path.abspath("."), relative_path)

# Dropdown menu to select between two images
def on_image_change(val):
    image_choice = val
    load_image(image_choice)
    ax_original.imshow(image)
    ax_saturation.imshow(image)
    ax_chroma.imshow(image)
    update_bar_graph()
    fig.canvas.draw_idle()

# List of images
images = ['andromedatry.png', 'RGB080604.png']

# Load default image
load_image(images[0])

# Initialize sample point for RGB bar graphs
sample_x, sample_y = image.shape[1] // 2, image.shape[0] // 2

# Create the main figure and axes
fig = plt.figure(figsize=(18, 12))
gs = fig.add_gridspec(3, 3, height_ratios=[3, 1, 1], hspace=0.3)

ax_original = fig.add_subplot(gs[0, 0])
ax_saturation = fig.add_subplot(gs[0, 1])
ax_chroma = fig.add_subplot(gs[0, 2])
bar_ax_original = fig.add_subplot(gs[1, 0])
bar_ax_saturation = fig.add_subplot(gs[1, 1])
bar_ax_chroma = fig.add_subplot(gs[1, 2])

# Display the original image
ax_original.imshow(image)
ax_original.set_title("Original")
ax_original.axis("off")

# Initialize saturation-adjusted and chroma-adjusted images
sat_image = image.copy()
chroma_image = image.copy()

# Display placeholders for saturation and chroma images
ax_saturation.imshow(sat_image)
ax_saturation.set_title("Saturation")
ax_saturation.axis("off")

ax_chroma.imshow(chroma_image)
ax_chroma.set_title("Chroma")
ax_chroma.axis("off")

# Initialize RGB bar graphs for each image
bar_ax_original.set_title("Original RGB Values")
bar_ax_original.set_ylim(0, 1)
sample_rgb = image[sample_y, sample_x].astype(np.float32) / 255.0  # Normalize to [0, 1]
bars_original = bar_ax_original.bar(['R', 'G', 'B'], sample_rgb, color=['red', 'green', 'blue'])

bar_ax_saturation.set_title("Saturation RGB Values")
bar_ax_saturation.set_ylim(0, 1)
bars_saturation = bar_ax_saturation.bar(['R', 'G', 'B'], sample_rgb, color=['red', 'green', 'blue'])

bar_ax_chroma.set_title("Chroma RGB Values")
bar_ax_chroma.set_ylim(0, 1)
bars_chroma = bar_ax_chroma.bar(['R', 'G', 'B'], sample_rgb, color=['red', 'green', 'blue'])

# Create a dropdown menu to select the image
from matplotlib.widgets import RadioButtons

radio_ax = plt.axes([0.35, 0.01, 0.3, 0.05])
radio = RadioButtons(radio_ax, images, active=0)
radio.on_clicked(on_image_change)

# Simplify slider appearance and callbacks
ax_sat_slider = plt.axes([0.35, 0.05, 0.2, 0.03], facecolor='lightgray')
sat_slider = Slider(ax_sat_slider, 'Saturation', 0.1, 2.0, valinit=1.0, valstep=0.01)
sat_slider.on_changed(on_slider_change)

ax_chroma_slider = plt.axes([0.7, 0.05, 0.2, 0.03], facecolor='lightgray')
chroma_slider = Slider(ax_chroma_slider, 'Chroma', 0.1, 2.0, valinit=1.0, valstep=0.01)
chroma_slider.on_changed(on_slider_change)

# Debounce timers
debounce_timer = None

# Connect sliders to debounce function
sat_slider.on_changed(on_slider_change)
chroma_slider.on_changed(on_slider_change)

# Connect mouse click to update sample point
fig.canvas.mpl_connect('button_release_event', on_slider_release)
fig.canvas.mpl_connect('button_press_event', onclick)

plt.show()
