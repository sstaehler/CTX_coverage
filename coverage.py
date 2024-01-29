#!/usr/bin/env python # -*- coding: utf-8 -*-
"""
Class for Spectral Fitting
:copyright:
    Simon Stähler (mail@simonstaehler.com), 2023
:license:
    None
"""
import numpy as np

RADIUS_MARS = 3389.5
LANDER_LAT = 4.5024
LANDER_LON = 135.6234


def area_annulus(radius_outer_km, radius_inner_km):
    theta_inner = radius_inner_km / RADIUS_MARS
    theta_outer = radius_outer_km / RADIUS_MARS

    area = 2 * np.pi * RADIUS_MARS**2 * (np.cos(theta_inner) - np.cos(theta_outer))
    return area


class Coverage:
    def __init__(self,
                 fnam_shape,
                 max_radius_km,
                 lat=LANDER_LAT,
                 lon=LANDER_LON):

        import geopandas as gpd
        import pandas as pd
        self.shapefile = gpd.read_file(fnam_shape)[
            ['CenterLat', 'CenterLon', 'MaxLat', 'MinLat', 'EastLon', 'WestLon', 'UTCstart',
             'UTCend', 'geometry', 'ProductId']]
        self.shapefile['date'] = pd.to_datetime(self.shapefile['UTCend'])
        self.lat = lat
        self.lon = lon
        circle = self.circle(radius_circle_km=max_radius_km)
        in_circle = circle.contains(self.shapefile["geometry"])
        self.shape_reduced = self.shapefile[in_circle]

    def in_annulus(self, radius_outer_km, radius_inner_km, geodataframe=None):
        if geodataframe is None:
            geodataframe = self.shape_reduced
        # in_distance = self.annulus(radius_outer_km, radius_inner_km).contains(geodataframe["geometry"])
        in_distance = geodataframe.intersection(other=self.annulus(radius_outer_km, radius_inner_km))
        return in_distance # geodataframe[in_distance]

    def area_images_in_annulus(self, radius_outer_km, radius_inner_km):
        joint = self.in_annulus(radius_outer_km, radius_inner_km).dissolve()
        return sum(joint.area * RADIUS_MARS)

    def area_images_in_annulus_time(self, radius_outer_km, radius_inner_km, time):
        images_in_time = self.shape_reduced.loc[self.shape_reduced["date"] > time].dissolve()
        # joint = self.in_annulus(radius_outer_km, radius_inner_km,
        #                         geodataframe=images_in_time).dissolve()
        joint = images_in_time.intersection(other=self.annulus(radius_outer_km, radius_inner_km))
        return sum(joint.area * RADIUS_MARS)

    def rel_imaged_area_in_annulus(self, radius_outer_km, radius_inner_km, time=None):
        if time is None:
            return self.area_images_in_annulus(radius_outer_km, radius_inner_km) / \
                area_annulus(radius_outer_km, radius_inner_km)
        else:
            return self.area_images_in_annulus_time(radius_outer_km, radius_inner_km, time) / \
                area_annulus(radius_outer_km, radius_inner_km)

    def circle(self, radius_circle_km):
        """
        Create polygon of circle with given radius
        :param radius_circle_km: float
        :return: tuple of polygon
        """
        import shapely.geometry
        poly = []
        for theta in range(0, 360, 2):
            poly.append(self.shoot(bearing_degree=theta, distance_km=radius_circle_km))
        polygon = shapely.geometry.Polygon(poly)
        return polygon

    def annulus(self, radius_outer_km, radius_inner_km):
        """
        Create polygon of annulus with given radius
        :param radius_outer_km: float
        :param radius_inner_km: float
        :return: tuple of polygon
        """
        import shapely.geometry
        poly_outer = []
        for theta in range(0, 361, 2):
            poly_outer.append(self.shoot(bearing_degree=theta, distance_km=radius_outer_km))
        poly_inner = []
        for theta in range(0, 361, 2):
            poly_inner.append(self.shoot(bearing_degree=theta, distance_km=radius_inner_km))
        polygon = shapely.geometry.Polygon(poly_outer, [poly_inner])
        return polygon

    def shoot(self, bearing_degree, distance_km, radius_km=3389.5):
        """
        Shoot a ray from point in direction for certain length and return where you land
        (Direct geodetic problem). Works on sphere
        :param latitude_1_degree: latitude of starting point
        :param longitude_1_degree: longitude of starting point
        :param bearing_degree: bearing from north, CW
        :param distance_km: distance in kilometer
        :param radius_km: radius of planet
        :return: latitude, longitude of target
        """
        lat1 = np.deg2rad(self.lat)
        lon1 = np.deg2rad(self.lon)
        bearing = np.deg2rad(bearing_degree)
        lat2 = np.arcsin(np.sin(lat1) * np.cos(distance_km / radius_km) +
                         np.cos(lat1) * np.sin(distance_km / radius_km) * np.cos(bearing))
        lon2 = lon1 + np.arctan2(np.sin(bearing) * np.sin(distance_km / radius_km) * np.cos(lat1),
                                 np.cos(distance_km / radius_km) - np.sin(lat1) * np.sin(lat2))
        # return np.rad2deg(lat2), np.mod(np.rad2deg(lon2) + 540., 360.) - 180.
        return np.rad2deg(lon2), np.rad2deg(lat2)
