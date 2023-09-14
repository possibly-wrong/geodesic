# Made with Natural Earth vector data @ naturalearthdata.com version 5.1.1.

import pygame
import numpy as np
import itertools
import os

water = (0, 0, 255, 255)
land = (0, 255, 0, 255)
border = (0, 0, 0, 255)

def xy(lonlat_deg):
    lon_deg, lat_deg = lonlat_deg
    return (int(width * (0.5 + lon_deg / 360)),
            int(height * (0.5 - lat_deg / 180)))

def draw_polygons(folder, image, color, fill=0):
    for filename in os.listdir(folder):
        lonlat = np.reshape(np.fromfile(os.path.join(folder, filename),
                                        dtype=np.float64), (-1, 2))
        polygons = [list(g)[:-1] for k, g in itertools.groupby(
            lonlat, key=lambda x: not np.isnan(x[0])) if k]
        for polygon in polygons:
            pygame.draw.polygon(image, color, [xy(p) for p in polygon], fill)

def make_map(filename, size, borders):
    global width, height
    width, height = size
    assert(width == 2 * height)
    image = pygame.Surface((width, height))
    image.fill(water)
    draw_polygons('ne_10m_land', image, land, 0)
    if borders:
        draw_polygons('ne_10m_admin_0_countries', image, border, 1)
        draw_polygons('ne_10m_admin_1_states_provinces', image, border, 1)
    pygame.image.save(image, filename)

if __name__ == '__main__':
    make_map('land_mask.png', (21600, 10800), False)
    make_map('world_map.png', (2700, 1350), True)
