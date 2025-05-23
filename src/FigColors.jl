using Colors

function fig_color()
    """
        fig_color()
    
    Outputs a palette of predefined colors for figures.
    """
    col = [
        RGB(30/255, 144/255, 205/255),      # Navy Blue
        RGB(255/255, 140/255, 0/255),       # Dark Orange
        RGB(50/255, 157/255, 13/255),       # Kinda Forest Green
        RGB(205/255, 175/255, 0/255),       # Gold yellow
        RGB(128/255, 0/255, 128/255),       # Purple
        RGB(205/255, 0/255, 0/255),         # Medium Red
        RGB(0/255, 0/255, 205/255),         # Medium Blue
        RGB(0/255, 205/255, 0/255),         # Medium green
        RGB(105/255, 105/255, 105/255),     # DimGrey1
        RGB(60/255, 60/255, 60/255)         # DimGrey2
    ]
    return col
end