# Example using geodesic.direct to "row east" along geodesics from east
# coast of North and South America.

import geodesic
import pygame
import math

water = (0, 0, 255, 255)
route = (255, 255, 0, 255)
lon_deg = -30

land_mask = pygame.image.load('land_mask.png')
world_map = pygame.image.load('world_map.png')
mask_width, mask_height = land_mask.get_size()
map_width, map_height = world_map.get_size()
image = pygame.Surface((map_width, map_height))
image.blit(world_map, (0, 0))

for lat_deg in range(-47, 64):
    x = int(mask_width * (0.5 + lon_deg / 360))
    y = int(mask_height * (0.5 - lat_deg / 180))
    while land_mask.get_at((x, y)) == water:
        x = x - 1
    x = x + 1
    lat1 = math.radians(90 - 180 * (y + 0.5) / mask_height)
    lon1 = math.radians(-180 + 360 * (x + 0.5) / mask_width)
    lat2, lon2 = lat1, lon1
    step = 0.5 * 2 * math.pi * 6378137 * math.cos(lat1) / mask_width
    steps = 0
    while land_mask.get_at((x, y)) == water:
        image.set_at((int(map_width * (0.5 + math.degrees(lon2) / 360)),
                      int(map_height * (0.5 - math.degrees(lat2) / 180))),
                     route)
        steps = steps + 1
        lat2, lon2 = geodesic.direct(lat1, lon1, math.pi / 2, step * steps)
        x = int(mask_width * (0.5 + math.degrees(lon2) / 360))
        y = int(mask_height * (0.5 - math.degrees(lat2) / 180))
pygame.image.save(image, 'east.png')
