import fitz  # PyMuPDF
from PIL import Image, ImageFont, ImageDraw, ImageChops
import os

# --- Define Project Root ---
projectRoot = './prj_190_viking_somalogic'

# --- Define Output Directory ---
outputDir = os.path.join(projectRoot, 'Scripts/Publication/supportingFiles/replication/locusZoom/output')

# Predefined gene names of interest
gene_names = ['AAMDC', 'PROCR', 'LTK', 'B3GAT1', 'BCL7A', 'COMMD10', 'NIF3L1', 'NTAQ1']  # Add more gene names here as needed

# Category mapping dictionary
categoryDict = {
    'B3GAT1_eQTLGen_rs78760579.pdf': ['base'],
    'B3GAT1_ieu-b-32_rs78760579.pdf': ['discovery', 'lymphocyte cell count'],
    'B3GAT1_ieu-b-85_rs78760579.pdf': ['discovery', 'prostate cancer'],
    'B3GAT1_ieu-b-4809_rs78760579.pdf': ['replication', 'prostate cancer'],
    'B3GAT1_Somalogic_rs78760579.pdf': ['base'],
    'B3GAT1_VIKING_rs78760579.pdf': ['base'],
    'LTK_ebi-a-GCST006867_rs1473781.pdf': ['discovery', 'type 2 diabetes'],
    'LTK_eQTLGen_rs1473781.pdf': ['base'],
    'LTK_Somalogic_rs1473781.pdf': ['base'],
    'LTK_type2_diabetes_rs1473781.pdf': ['replication', 'type 2 diabetes'],
    'LTK_ukb-b-10753_rs1473781.pdf': ['discovery', 'type 2 diabetes'],
    'LTK_ukb-b-14609_rs1473781.pdf': ['discovery', 'type 2 diabetes'],
    'LTK_VIKING_rs1473781.pdf': ['base'],
    'NIF3L1_ebi-a-GCST010723_rs10931931.pdf': ['discovery', 'early age-related macular degeneration'],
    'NIF3L1_eQTLGen_rs10931931.pdf': ['base'],
    'NIF3L1_Somalogic_rs10931931.pdf': ['base'],
    'NIF3L1_VIKING_rs10931931.pdf': ['base'],
    'NTAQ1_eQTLGen_rs13258747.pdf': ['base'],
    'NTAQ1_Somalogic_rs13258747.pdf': ['base'],
    'NTAQ1_VIKING_rs13258747.pdf': ['base'],
    'NTAQ1_ieu-b-4865_rs13258747.pdf': ['discovery', 'testosterone levels'],
    'AAMDC_ebi-a-GCST004599_rs72941336.pdf': ['discovery', 'mean platelet volume'],
    'AAMDC_eQTLGen_rs72941336.pdf': ['base'],
    'AAMDC_GCST90014290_rs72941336.pdf': ['discovery', 'intrinsic epigenetic age acceleration'],
    'AAMDC_GCST90244093_rs72941336.pdf': ['replication', 'forced vital capacity'],
    'AAMDC_GCST90310295_rs72941336.pdf': ['replication', 'diastolic blood pressure'],
    'AAMDC_ieu-a-1006_rs72941336.pdf': ['replication', 'mean platelet volume'],
    'AAMDC_ieu-b-39_rs72941336.pdf': ['discovery', 'diastolic blood pressure'],
    'AAMDC_ieu-b-105_rs72941336.pdf': ['discovery', 'forced vital capacity'],
    'AAMDC_Olink_rs72941336.pdf': ['base'],
    'AAMDC_Somalogic_rs72941336.pdf': ['base'],
    'AAMDC_ukb-b-7953_rs72941336.pdf': ['discovery', 'forced vital capacity'],
    'AAMDC_ukb-b-7992_rs72941336.pdf': ['discovery', 'diastolic blood pressure'],
    'AAMDC_VIKING_rs72941336.pdf': ['base'],
    'BCL7A_eQTLGen_rs1169084.pdf': ['base'],
    'BCL7A_GCST90310294_rs1169084.pdf': ['replication', 'systolic blood pressure'],
    'BCL7A_ieu-b-38_rs1169084.pdf': ['discovery', 'systolic blood pressure'],
    'BCL7A_Olink_rs1169084.pdf': ['base'],
    'BCL7A_Somalogic_rs1169084.pdf': ['base'],
    'BCL7A_VIKING_rs1169084.pdf': ['base'],
    'COMMD10_ebi-a-GCST006696_rs56953556.pdf': ['discovery', 'parental longevity (mother\'s age)'],
    'COMMD10_eQTLGen_rs56953556.pdf': ['base'],
    'COMMD10_maternal_longevity_eLife_Timmers_rs56953556.pdf': ['replication', 'parental longevity (mother\'s age)'],
    'COMMD10_Somalogic_rs56953556.pdf': ['base'],
    'COMMD10_VIKING_rs56953556.pdf': ['base'],
    'PROCR_BMI_GIANT_rs6060241.pdf': ['replication', 'body mass index'],
    'PROCR_eQTLGen_rs6060241.pdf': ['base'],
    'PROCR_Olink_rs6060241.pdf': ['base'],
    'PROCR_Somalogic_rs6060241.pdf': ['base'],
    'PROCR_ukb-b-8909_rs6060241.pdf': ['discovery', 'body mass index'],
    'PROCR_ukb-b-12854_rs6060241.pdf': ['discovery', 'body mass index'],
    'PROCR_ukb-b-19953_rs6060241.pdf': ['discovery', 'body mass index'],
    'PROCR_ukb-b-20531_rs6060241.pdf': ['discovery', 'body mass index'],
    'PROCR_VIKING_rs6060241.pdf': ['base']
    # Add more filenames and their corresponding categories
}
base_subcategories = ['VIKING', 'Somalogic', 'Olink', 'eQTLGen']

def extract_gene_name(filename):
    """Extract the gene name from the filename."""
    for gene in gene_names:
        if gene in filename:
            return gene
    return None  # Return None if no gene is found

def extract_base_subcategory(filename):
    """Extract the base subcategory from the filename."""
    for subcategory in base_subcategories:
        if subcategory in filename:
            return subcategory
    return None  # Return None if no subcategory is found

def create_label_image(label, width, height=100, bg_color=(200, 200, 200)):
    """Create an image with a label (text) centered horizontally and a customizable background color."""
    label_img = Image.new('RGB', (width, height), bg_color)  # Background color is passed as an argument
    draw = ImageDraw.Draw(label_img)

    # Try to load the DejaVuSans font, which is commonly available and can be resized
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 70)
    except IOError:
        # If DejaVuSans is not found, provide a fallback to a smaller, built-in font (but not scalable)
        print("DejaVuSans font not found. Using default small font, which may be too small.")
        font = ImageFont.load_default()

    # Calculate text size and position
    text_bbox = draw.textbbox((0, 0), label, font=font)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]
    text_x = (width - text_width) // 2  # Center horizontally
    text_y = (height - text_height) // 2  # Center vertically
    
    # Draw the label text
    draw.text((text_x, text_y), label, font=font, fill=(0, 0, 0))  # Black text
    
    return label_img

def overlay_color_on_image(image, color, alpha=0.4):
    """Overlay a transparent color on top of an image."""
    overlay = Image.new('RGB', image.size, color)
    return Image.blend(image, overlay, alpha)


def render_first_page_as_image(pdf_file, dpi=200, title=None):
    """Render the first page of a PDF as an image and overlay a title at the top right."""
    pdf_path = os.path.join(outputDir, pdf_file)
    
    try:
        doc = fitz.open(pdf_path)
        page = doc.load_page(0)  # Get the first page
        
        # Render page to an image with the desired resolution (DPI)
        zoom = dpi / 72  # PyMuPDF uses 72dpi as the base
        mat = fitz.Matrix(zoom, zoom)  # Scale the page by the zoom factor
        
        # Render the page as a pixmap (image)
        pix = page.get_pixmap(matrix=mat)
        
        # Convert the pixmap to a PIL Image
        img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
        
        # Load a font for the title
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 50)
        
        # If a title is provided, overlay it
        if not title:
            title = ''

        #Figure out why all titles are overriden by none
        #Figure out why some filenames are not plotted like COMMD10_maternal_longevity_eLife_Timmers_rs56953556.pdf
        draw = ImageDraw.Draw(img)
        draw.text((400, 80), title, fill=(255, 0, 0), font = font)

        print(f"Page rendered as image from {pdf_file} with title: {title}")
        return img
    
    except Exception as e:
        print(f"Error processing {pdf_file}: {e}")
    
    return None



def create_vertical_label_image(labels, height_per_label=100, width=300, bg_color=(200, 200, 200)):
    """Create a vertical label image to display labels on the left side of the image with specified width."""
    total_height = len(labels) * height_per_label
    label_img = Image.new('RGB', (width, total_height), bg_color)  # Background color for labels
    draw = ImageDraw.Draw(label_img)

    try:
        # Try to load a font that can be resized
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 48)  # Adjusted font size
    except IOError:
        font = ImageFont.load_default()

    # Draw each label vertically aligned in the label image
    for i, label in enumerate(labels):
        # Calculate the width of the text
        text_bbox = draw.textbbox((0, 0), label, font=font)
        text_width = text_bbox[2] - text_bbox[0]
        text_height = text_bbox[3] - text_bbox[1]

        # If text width exceeds the available width, reduce font size until it fits
        while text_width > width:
            font = ImageFont.truetype("DejaVuSans-Bold.ttf", font.size - 1)  # Reduce font size
            text_bbox = draw.textbbox((0, 0), label, font=font)
            text_width = text_bbox[2] - text_bbox[0]
            text_height = text_bbox[3] - text_bbox[1]

        # Recalculate position for centering the text horizontally
        text_x = (width - text_width) // 2  # Center horizontally
        text_y = (height_per_label - text_height) // 2 + i * height_per_label  # Center vertically

        # Draw the text
        draw.text((text_x, text_y), label, font=font, fill=(0, 0, 0))  # Black text
    
    return label_img


def create_placeholder_image(width, height):
    """Create a blank placeholder image for missing subcategories."""
    return Image.new('RGB', (width, height), (255, 255, 255))  # White blank image

def process_base_images(image_dict, gene, placeholder_width, placeholder_height):
    """Process the base subcategory images for a given gene."""
    base_images = []
    for subcategory in base_subcategories:
        if gene in image_dict and 'base' in image_dict[gene] and subcategory in image_dict[gene]['base']:
            img = image_dict[gene]['base'][subcategory]
            base_images.append(img)
        else:
            # Create a placeholder that says "Not measured" for missing images
            base_images.append(create_placeholder_image(placeholder_width, placeholder_height, text="Not measured"))
    return base_images

def process_subgroups(gene, image_dict):
    """Process the discovery and replication subgroups for a given gene."""
    subgroups = {}
    for pdf_file, categories in categoryDict.items():
        if gene in pdf_file:
            if len(categories) > 1:
                subgroup_name = categories[1]
            else:
                continue  # Skip if no valid subgroup name

            if subgroup_name not in subgroups:
                subgroups[subgroup_name] = {'discovery': [], 'replication': []}

            #Get study name from pdf_file
            title = pdf_file.split('_')[1:-1]
            title = '_'.join(title)

            if 'discovery' in categories:
                img = render_first_page_as_image(pdf_file, title=title)
                if img:
                    subgroups[subgroup_name]['discovery'].append(img)
            elif 'replication' in categories:
                img = render_first_page_as_image(pdf_file, title=title)
                if img:
                    subgroups[subgroup_name]['replication'].append(img)

    return subgroups

def stack_subgroup_images(subgroups, placeholder_width):
    """Stack images within each subgroup: discovery above replication."""
    subgroup_stacked_images = []
    for subgroup, images in subgroups.items():
        subgroup_images = []

        # Add a label for the subgroup (e.g., 'lymphocyte cell count')
        label_image = create_label_image(subgroup, placeholder_width, 100, (200, 200, 200))
        subgroup_images.append(label_image)

        # Add discovery images with green overlay
        if images['discovery']:
            for img in images['discovery']:
                if img:
                    img_with_overlay = overlay_color_on_image(img, (180, 238, 180), alpha=0.2)  # Softer light green
                    subgroup_images.append(img_with_overlay)

        # Add replication images with blue overlay
        if images['replication']:
            for img in images['replication']:
                if img:
                    img_with_overlay = overlay_color_on_image(img, (173, 216, 230), alpha=0.4)  # Light blue
                    subgroup_images.append(img_with_overlay)

        # Stack discovery, replication, and label images within the subgroup
        if subgroup_images:
            total_height = sum(img.height for img in subgroup_images)
            stacked_subgroup_image = Image.new('RGB', (subgroup_images[0].width, total_height))
            y_offset = 0
            for img in subgroup_images:
                stacked_subgroup_image.paste(img, (0, y_offset))
                y_offset += img.height
            subgroup_stacked_images.append(stacked_subgroup_image)

    return subgroup_stacked_images

def stack_images_vertically(image_dict, output_file, separator_thickness=30):
    columns = []

    for gene in gene_names:
        gene_images = []

        placeholder_width, placeholder_height = get_placeholder_dimensions(image_dict, gene)

        base_images = process_base_images(image_dict, gene, placeholder_width, placeholder_height)
        if base_images:
            stacked_base_image = create_hatched_background(
                base_images[0].width,
                sum(img.height for img in base_images)
            )
            y_offset = 0
            for img in base_images:
                stacked_base_image.paste(img, (0, y_offset))
                y_offset += img.height
            gene_images.append(stacked_base_image)

        subgroups = process_subgroups(gene, image_dict)
        subgroup_images = stack_subgroup_images(subgroups, placeholder_width)
        gene_images.extend(subgroup_images)

        if gene_images:
            gene_label_image = create_label_image(gene, placeholder_width, 100)
            gene_images.insert(0, gene_label_image)

            gene_column = stack_gene_images_with_separator(gene_images, separator_thickness)
            columns.append(gene_column)

    if not columns:
        raise RuntimeError("No columns generated")

    # Using the placeholder dimensions as a guide
    legend_w = int(placeholder_width * 1)
    legend_h = int(placeholder_height * 1)
    legend_image = create_legend_image(width=legend_w, height=legend_h)

    # Make the legend the last frame of the last column (bottom-right of canvas)
    columns = append_legend_to_last_column(columns, legend_image, separator_thickness)

    # Build the left-side labels after we know placeholder_height
    base_labels = base_subcategories
    label_image = create_vertical_label_image(base_labels, width=200, height_per_label=placeholder_height)

    # Stack columns + labels into the final image
    final_image = stack_columns_with_labels(label_image, columns)

    final_image = add_global_diagonal_watermark(
        final_image,
        label_width=label_image.width,
        label_height=label_image.height,
        columns=columns,
        text="Empty space\n(no measurements)",
        angle=35,
        opacity=80
    )

    final_image.save(output_file)
    print(f"Image saved to {output_file}")


def get_placeholder_dimensions(image_dict, gene):
    """Get the placeholder dimensions for a given gene."""
    for subcategory in base_subcategories:
        if gene in image_dict and 'base' in image_dict[gene] and subcategory in image_dict[gene]['base']:
            img = image_dict[gene]['base'][subcategory]
            return img.width, img.height
    return 500, 500  # Default size if no image is found

def create_hatched_background(width, height, bg_color=(220, 220, 220), line_color=(200, 200, 200), spacing=20):
    img = Image.new('RGB', (width, height), bg_color)
    draw = ImageDraw.Draw(img)
    for x in range(-height, width, spacing):
        draw.line((x, 0, x + height, height), fill=line_color)
    return img

def stack_gene_images_with_separator(gene_images, separator_thickness):

    """Stack gene images with separators."""
    width = gene_images[0].width
    separator = create_hatched_background(width, separator_thickness)  # Light grey hatched separator

    total_height = sum(img.height for img in gene_images) + separator_thickness * (len(gene_images) - 1)
    stacked_gene_image = create_hatched_background(width, total_height)
    
    y_offset = 0
    for i, img in enumerate(gene_images):
        stacked_gene_image.paste(img, (0, y_offset))
        y_offset += img.height
        if i < len(gene_images) - 1:
            stacked_gene_image.paste(separator, (0, y_offset))
            y_offset += separator_thickness

    return stacked_gene_image

def stack_columns_with_labels(label_image, columns):
    total_width = sum(col.width for col in columns) + label_image.width
    max_height = max(col.height for col in columns)
    total_height = max(max_height, label_image.height)

    final_image = create_hatched_background(total_width, total_height)

    # Paste label image on the left
    final_image.paste(label_image, (0, 0))

    # Paste columns to the right of labels
    x_offset = label_image.width
    for col in columns:
        final_image.paste(col, (x_offset, 0))
        x_offset += col.width

    return final_image

def create_placeholder_image(width, height, text="Not measured"):
    """Create an image with a 'Not measured' label when an image is missing."""
    placeholder_img = create_hatched_background(width, height)
    draw = ImageDraw.Draw(placeholder_img)

    try:
        # Load a font that can be resized
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 48)  # Adjusted font size
    except IOError:
        font = ImageFont.load_default()

    # Calculate text size and position
    text_bbox = draw.textbbox((0, 0), text, font=font)
    text_width = text_bbox[2] - text_bbox[0]
    text_height = text_bbox[3] - text_bbox[1]
    text_x = (width - text_width) // 2  # Center horizontally
    text_y = (height - text_height) // 2  # Center vertically

    # Draw the text "Not measured"
    draw.text((text_x, text_y), text, font=font, fill=(0, 0, 0))  # Black text
    
    return placeholder_img

def create_legend_image(width=400, height=150, bg_color=(255, 255, 255)):
    """Create a legend image that explains the color palette used with color boxes enclosed in a black outline."""
    legend_img = Image.new('RGB', (width, height), bg_color)  # White background
    draw = ImageDraw.Draw(legend_img)

    # Define the text labels and colors (without the color names)
    labels = [
        ("molecular measurement GWAS", (255, 255, 255)),  # White background
        ("MR discovery GWAS", (180, 238, 180)),           # Green for discovery GWAS
        ("MR replication GWAS", (173, 216, 230))          # Blue for replication GWAS
    ]

    # Load font
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", 80)  # Use a larger font for the legend
    except IOError:
        font = ImageFont.load_default()

    # Define positions and rectangle sizes
    padding = 50
    box_size = 200  # Size of the color box
    text_x_offset = box_size + padding * 2  # Increase offset to add more space between box and text
    y_offset = padding  # Start offset from the top

    outline_thickness = 5  # Thickness of the black outline

    for description, color in labels:
        # Draw the black outline first (slightly larger than the color box)
        draw.rectangle(
            [padding - outline_thickness, y_offset - outline_thickness,
             padding + box_size + outline_thickness, y_offset + box_size + outline_thickness],
            outline="black", width=outline_thickness
        )

        # Draw the color box
        draw.rectangle([padding, y_offset, padding + box_size, y_offset + box_size], fill=color)

        # Draw the description text (without the color names)
        draw.text((text_x_offset, y_offset), description, font=font, fill=(0, 0, 0))

        # Move the y_offset down for the next label
        y_offset += box_size + padding

    return legend_img

def append_legend_to_last_column(columns, legend_image, separator_thickness):
    """
    Make the legend the last 'frame' of the last column.
    Pad the last column so the legend sits at the canvas bottom,
    and right-align the legend within the column.
    """
    if not columns:
        return columns

    last_col = columns[-1]
    tallest = max(c.height for c in columns)

    needed_height = tallest - last_col.height - separator_thickness - legend_image.height
    spacer_h = max(0, needed_height)

    new_h = last_col.height + spacer_h + separator_thickness + legend_image.height
    new_w = last_col.width

    # Use HATCHED background for the rebuilt last column
    new_last = create_hatched_background(new_w, new_h)
    y = 0
    new_last.paste(last_col, (0, y)); y += last_col.height

    if spacer_h > 0:
        # Use HATCHED spacer too
        spacer = create_hatched_background(new_w, spacer_h)
        new_last.paste(spacer, (0, y)); y += spacer_h

    # Hatched separator already via helper
    sep = create_hatched_background(new_w, separator_thickness)
    new_last.paste(sep, (0, y)); y += separator_thickness

    # right-align legend within the column
    legend_x = new_w - legend_image.width
    new_last.paste(legend_image, (legend_x, y))

    columns[-1] = new_last
    return columns

def add_global_diagonal_watermark(final_image, label_width, label_height, columns, text,
                                  angle=35, opacity=80):
    """
    Overlay one large diagonal watermark across *only* the empty grey/hatched
    areas (all bottom rows under shorter columns and under the label panel).
    - final_image: PIL.Image (RGB)
    - label_width/label_height: size of the left label area you pasted
    - columns: list of column images in the order you pasted them
    - text: watermark string
    - angle: rotation in degrees
    - opacity: 0..255 (per-pixel max opacity for watermark)
    """
    W, H = final_image.size
    tallest = max(c.height for c in columns)

    # 1) Build an L-mode mask marking "empty areas" as 255, else 0
    empty_mask = Image.new('L', (W, H), 0)
    mdraw = ImageDraw.Draw(empty_mask)

    # Empty under the left label block
    if label_height < tallest:
        mdraw.rectangle([0, label_height, label_width, tallest], fill=255)

    # Empty under each column
    x = label_width
    for col in columns:
        if col.height < tallest:
            mdraw.rectangle([x, col.height, x + col.width, tallest], fill=255)
        x += col.width

    # 2) Create a big text layer (same size as final image)
    text_layer = Image.new('RGBA', (W, H), (0, 0, 0, 0))
    draw = ImageDraw.Draw(text_layer)

    # Choose a font size roughly proportional to diagonal length
    diag = (W**2 + tallest**2) ** 0.5
    target_size = max(14, int(diag * 0.04))
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", target_size)
    except IOError:
        font = ImageFont.load_default()

    # Draw text centered before rotation
    tb = draw.textbbox((0, 0), text, font=font)
    tw, th = tb[2] - tb[0], tb[3] - tb[1]
    # Shift factors: fraction of available space from top-left
    shift_x = 0.45  # 0.5 would be centered horizontally
    shift_y = 0.65  # 0.5 would be centered vertically

    cx = int(W * shift_x)
    cy = int(tallest * shift_y)
    tx = cx - tw // 2
    ty = cy - th // 2
    draw.text((tx, ty), text, font=font, fill=(0, 0, 0, opacity))

    # Rotate around the canvas center without changing size
    text_layer = text_layer.rotate(angle, expand=False, resample=Image.BICUBIC, center=(W/2, H/2))

    # 3) Restrict watermark to the "empty areas" and keep its original glyph alpha
    glyph_alpha = text_layer.split()[3]  # A channel of glyphs
    # Scale the empty mask to (0..opacity) and multiply with glyph alpha
    scaled_mask = empty_mask.point(lambda p: int(p * opacity / 255))
    combined_alpha = ImageChops.multiply(glyph_alpha, scaled_mask)
    text_layer.putalpha(combined_alpha)

    # 4) Composite and return
    out = Image.alpha_composite(final_image.convert('RGBA'), text_layer).convert('RGB')
    return out



# Dictionary to hold images sorted by gene and category, with subcategories in base
image_dict = {gene: {'base': {}} for gene in gene_names}

# Process each PDF
fileNames = [f for f in os.listdir(outputDir) if f.endswith('.pdf')]

#Process base images (VIKING, Somalogic, Olink, eQTLGen)
for pdf_file in fileNames:
    gene = extract_gene_name(pdf_file)
    if gene:
        category = categoryDict.get(pdf_file, None)
        if 'base' in category:
            subcategory = extract_base_subcategory(pdf_file)
            if subcategory:
                img = render_first_page_as_image(pdf_file)
                if img:
                    # Save the image in the dictionary under the correct gene and subcategory
                    image_dict[gene]['base'][subcategory] = img

# Now process the discovery/replication & stack the images vertically per gene with a separator and save the final image
outputImageFile = os.path.join(outputDir, 'stacked_images.png')
stack_images_vertically(image_dict, outputImageFile)