"""Unitary tests for the eponymous script.

Copyright (c) 2024 Regional-Modeling-LATMOS-IGE.

This software is released under the terms of the Expat (aka MIT) license:

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

from pytest import raises
import numpy as np
try:
    from shapely import Polygon
except ImportError:
    from shapely.geometry import Polygon
import cams2wrfchem as tt

## unique_index

a = np.array([4, 7, 21, 5, -8, 44, 2, -5, 8, 21, 10, -100])
b = np.array([[4, 7, 21, 5], [-8, 44, 2, -5], [8, 21, 6, 2]])

def test__unique_index__expected_use_cases():
    assert tt.unique_index(a, 4) == 0
    assert tt.unique_index(a, 5) == 3
    assert tt.unique_index(a, 44) == 5

def test__unique_index__bad_array_shape():
    with raises(ValueError):
        tt.unique_index(b, 4)

def test__unique_index__inexisting_value():
    with raises(ValueError):
        tt.unique_index(a, 999)

def test__unique_index__repeated_value():
    with raises(ValueError):
        tt.unique_index(a, 21)

## calc_cellweights (optimized and brute force versions)

x_bdy = np.array([list(range(6)) for i in range(7)])
y_bdy = np.array([[i]*6 for i in range(6, -1, -1)])
x = x_bdy[:-1,:-1] + 0.5
y = y_bdy[:-1,:-1] - 0.5

def test__calc_cellweights__regular_square_grid():
    # Be careful, the order of the elements in the answer should not matter
    polygon = Polygon([(1,2), (1,3), (2,3), (2,2)])
    answer = [tt.CellWeight(3, 1, 1)]
    assert tt.calc_cellweights(x, y, x_bdy, y_bdy, polygon) == answer
    polygon = Polygon([(1,2), (1,3), (3,3), (3,2)])
    answer = [tt.CellWeight(3, 1, 0.5), tt.CellWeight(3, 2, 0.5)]
    assert set(tt.calc_cellweights(x, y, x_bdy, y_bdy, polygon)) == set(answer)
    polygon = Polygon([(0,4), (0,6), (2.5,6), (2.5,4)])
    answer = [tt.CellWeight(0, 0, 0.2), tt.CellWeight(1, 0, 0.2),
              tt.CellWeight(0, 1, 0.2), tt.CellWeight(1, 1, 0.2),
              tt.CellWeight(0, 2, 0.1), tt.CellWeight(1, 2, 0.1)]
    assert set(tt.calc_cellweights(x, y, x_bdy, y_bdy, polygon)) == set(answer)

def test__calc_cellweights_bruteforce__regular_square_grid():
    # Be careful, the order of the elements in the answer should not matter
    polygon = Polygon([(1,2), (1,3), (2,3), (2,2)])
    function = tt.calc_cellweights_bruteforce
    answer = [tt.CellWeight(3, 1, 1)]
    assert function(x_bdy, y_bdy, polygon) == answer
    polygon = Polygon([(1,2), (1,3), (3,3), (3,2)])
    answer = [tt.CellWeight(3, 1, 0.5), tt.CellWeight(3, 2, 0.5)]
    assert set(function(x_bdy, y_bdy, polygon)) == set(answer)
    polygon = Polygon([(0,4), (0,6), (2.5,6), (2.5,4)])
    answer = [tt.CellWeight(0, 0, 0.2), tt.CellWeight(1, 0, 0.2),
              tt.CellWeight(0, 1, 0.2), tt.CellWeight(1, 1, 0.2),
              tt.CellWeight(0, 2, 0.1), tt.CellWeight(1, 2, 0.1)]
    assert set(function(x_bdy, y_bdy, polygon)) == set(answer)
