import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data_frame_names = [
    "step",
    "time",
    "area",
    "area-x",
    "area-y",
    "area-z",
    "f-conv-x",
    "f-conv-y",
    "f-conv-z",
    "f-p-x",
    "f-p-y",
    "f-p-z",
    "f-visc-x",
    "f-visc-y",
    "f-visc-z",
    "f-body-x",
    "f-body-y",
    "f-body-z",
    "f-total-x",
    "f-total-y",
    "f-total-z"
]

pathname = r"C:\Users\gusro\OneDrive\Desktop\ME123\ME123 Final Project\Default Windmills\data\force_islands.dat"

base_island = pd.read_csv(pathname,
                            delimiter='\s+', engine='python', skiprows=3, encoding='utf-8',names=data_frame_names)