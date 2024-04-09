#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 14:34:30 2024

@author: raharinirina
"""

def Generate_txt(list_items, file):
    tofile = open(file, "w")
    for i in range(len(list_items)):
        tofile.write(list_items[i]+"\n")
    tofile.close()
        
    