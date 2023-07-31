import json
import colorsys
from matplotlib import colors as mc
import operator

file = open("colorbrewer.json","r")
colorbrewer = json.load(file)
file.close()

colormaps = [u'Oranges',u'Blues',u'Reds',u'Purples',u'Greens',u'Greys',u'Oranges',u'Blues',u'Reds',u'Purples',u'Greens',u'Greys',u'Oranges',u'Blues',u'Reds',u'Purples',u'Greens',u'Greys']

colormaps = [u'Blues',u'Oranges',u'Greens',u'Reds',u'Purples',u'Greys',u'Blues',u'Oranges',u'Greens',u'Reds',u'Purples',u'Greys',u'Blues',u'Oranges',u'Greens',u'Reds',u'Purples',u'Greys']

colormaps = [u'Blues',u'Oranges',u'Greens',u'Purples',u'Reds',u'Blues',u'Oranges',u'Greens',u'Purples',u'Reds',u'Blues',u'Oranges',u'Greens',u'Purples',u'Reds']


def parse_colorstring(colorstring):
    
    items = colorstring[4:-1].split(",")
    
    hex_items = []
    for item in items:
        hex_str = hex(int(item))[2:]
        #print hex_str
        if len(hex_str)<2:
            hex_str = '0'+hex_str
        hex_items.append(hex_str)
    
    color = ('#'+"".join(hex_items)).upper()
    return color
    
def get_colors(colormap_idx,num):
    
    double=False
    if num<3:
        num=3
    if num>9:
        num=9
        double=True
        
    colors = [parse_colorstring(item) for item in colorbrewer[colormaps[colormap_idx]][str(num)]]
    
    color_items = [(color,colorsys.rgb_to_hls(*mc.to_rgb(color))[1]) for color in colors]
    
    
    sorted_color_items = list(sorted(color_items, key=operator.itemgetter(1),reverse=False))
     
    colors = [item[0] for item in sorted_color_items]
    
    if double:
        colors = colors+colors
    return colors
    
if __name__=='__main__':
    print(colorbrewer[u'Blues']['4'])