using Colors, ColorSchemes, CairoMakie

function update_fig_theme()
    fs = 14 # FontSize
    fig_theme = Theme(
        #palette = (color = fig_color(), ),
        fontsize = fs,
        figure_padding = 8,
        rowgap = 8,
        colgap = 8,
        Axis = (xlabelfont = "Times", ylabelfont = "Times",
            xticklabelfont = "Times", yticklabelfont = "Times",
            xgridvisible = false, ygridvisible = false,
            xtickalign = 1, ytickalign = 1,
            xminorticksvisible = true, yminorticksvisible = true,
            xminortickalign = 1, yminortickalign = 1),
        Legend = (framecolor = :transparent, backgroundcolor = :transparent,
            patchsize = (30, 10),
            labelfont = "Times", labelsize = fs),
        Text = (font = "Times", fontsize = fs),
        palette = (color = ColorSchemes.seaborn_muted, marker = ColorSchemes.seaborn_muted),
        #palette = (color = fig_color(), marker = fig_color())
        #https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/
    )
    set_theme!(fig_theme)
    return nothing
end

function fig_color()
    """
        fig_color()

    Outputs a palette of predefined colors for figures. seaborn_muted or tableau_10
    #https://juliagraphics.github.io/ColorSchemes.jl/stable/catalogue/
    Idenfity the color with:
    using Colors, ColorSchemes
    # One color
    rgb_255 = round.(Int, 255 .* collect((red(col[1]), green(col[1]), blue(col[1]))))
    # All colors
    rgb_255 = Matrix{Integer}(undef, 3, length(col))
    for (i, j) in enumerate(col)
        rgb_255[:, i] = round.(Int, 255 .* collect((red(j), green(j), blue(j))))
    end
    # Colors for Seaborn Color Blind:
    1	115	178     # Blue
    222	143	5       # Yellow
    2	158	115     # Green
    213	94	0       # Red
    204	120	188     # Purple
    202	145	97      # Brown
    251	175	228     # Pink
    148	148	148     # Grey
    236	225	51      # Yellow (light)
    86	180	233     # Blue (light)
    """
    col = ColorSchemes.seaborn_colorblind
    #col = ColorSchemes.tableau_10
    return col
end

function fig_color1()
    """
        fig_color()

    Outputs a palette of predefined colors for figures. Fresh look with muted colors.
    See: https://projects.susielu.com/viz-palette
    https://www.figma.com/color-wheel/
    """
    col = [
        RGB(40 / 255, 100 / 255, 180 / 255),        # 1 - Blue (Complemetary 1)
        RGB(220 / 255, 120 / 255, 20 / 255),        # 2 - Orange (Complemetary 1)
        RGB(190 / 255, 150 / 255, 0 / 255),         # 3 - Yellow (Complemetary 2)
        RGB(180 / 255, 40 / 255, 110 / 255),        # 4 - Purple (Complemetary 2)
        RGB(100 / 255, 180 / 255, 20 / 255),        # 5 - Green (Complemetary 3)
        RGB(180 / 255, 30 / 255, 30 / 255),         # 6 - Red (Complemetary 3)
        RGB(150 / 255, 150 / 255, 150 / 255),       # 7 - Light Grey
        RGB(100 / 255, 100 / 255, 100 / 255),       # 8 - Dim Grey
        RGB(50 / 255, 50 / 255, 50 / 255)           # 9 - Dark Grey
    ]
    return col
end

function fig_color2()
    """
        fig_color2()

    Outputs a palette of predefined colors for figures. Brighter for previous articles.
    """
    col = [
        RGB(30 / 255, 144 / 255, 205 / 255),      # 1 - Navy Blue
        RGB(255 / 255, 140 / 255, 0 / 255),       # 2 - Dark Orange
        RGB(50 / 255, 157 / 255, 13 / 255),       # 3 - Kinda Forest Green
        RGB(205 / 255, 175 / 255, 0 / 255),       # 4 - Gold yellow
        RGB(168 / 255, 40 / 255, 168 / 255),      # 5 - Purple
        RGB(205 / 255, 0 / 255, 0 / 255),         # 6 - Medium Red
        RGB(0 / 255, 0 / 255, 205 / 255),         # 7 - Medium Blue
        RGB(0 / 255, 205 / 255, 0 / 255),         # 8 - Medium green
        RGB(105 / 255, 105 / 255, 105 / 255),     # 9 - DimGrey1
        RGB(60 / 255, 60 / 255, 60 / 255)         # 10 - DimGrey2
    ]
    return col
end

function fig_color3()
    """
        fig_color()

    Outputs a palette of predefined colors for figures. Bit darker.
    """
    col = [
        RGB(10 / 255, 120 / 255, 190 / 255),      # 1 - Navy Blue
        RGB(240 / 255, 130 / 255, 0 / 255),       # 2 - Dark Orange
        RGB(30 / 255, 150 / 255, 10 / 255),       # 3 - Kinda Forest Green
        RGB(190 / 255, 150 / 255, 0 / 255),       # 4 - Gold yellow
        RGB(150 / 255, 20 / 255, 150 / 255),      # 5 - Purple
        RGB(180 / 255, 20 / 255, 20 / 255),       # 6 - Medium Red
        RGB(20 / 255, 20 / 255, 180 / 255),       # 7 - Medium Blue
        RGB(20 / 255, 180 / 255, 20 / 255),       # 8 - Medium green
        RGB(150 / 255, 150 / 255, 150 / 255),     # 9 - DimGrey1
        RGB(100 / 255, 100 / 255, 100 / 255),     # 10 - DimGrey2
        RGB(50 / 255, 50 / 255, 50 / 255)         # 10 - DimGrey2
    ]
    return col
end

function fig_color4()
    """
        fig_color()

    Outputs a palette of predefined colors for figures.
    Matlab "gem" 2025a colors: https://www.mathworks.com/help/matlab/ref/colororder.html
    """
    col = [
        RGB(0.0660, 0.4430, 0.7450),            # 1 - Blue (Complemetary 1)
        RGB(0.8660, 0.3290, 0.0000),            # 2 - Orange (Complemetary 1)
        RGB(0.9290, 0.6940, 0.1250),            # 3 - Yellow (Complemetary 2)
        RGB(0.5210, 0.0860, 0.8190),            # 4 - Purple (Complemetary 2)
        RGB(0.2310, 0.6660, 0.1960),            # 5 - Green (Complemetary 3)
        RGB(0.8190, 0.0150, 0.5450),            # 6 - Red (Complemetary 3)
        RGB(150 / 255, 150 / 255, 150 / 255),       # 7 - Light Grey
        RGB(100 / 255, 100 / 255, 100 / 255),       # 8 - Dim Grey
        RGB(50 / 255, 50 / 255, 50 / 255)           # 9 - Dark Grey
    ]
    return col
end