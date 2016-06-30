;Extrapolate the position of a point along a line surrounding it
FUNCTION fit_lin, x_arr, y_arr, x
temp = max(where(x_arr LT x),i0)
i1 = i0 + 1
slope = (y_arr[i1] - y_arr[i0])/(x_arr[i1] - x_arr[i0])
y = slope*(x - x_arr[i0]) + y_arr[i0]
RETURN, y
END
