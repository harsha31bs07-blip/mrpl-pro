* Small LP used by the reader tests.
*   max  3 x + 2 y
*   s.t. x +   y <= 4
*        x + 3 y <= 6
*        x       <= 3
*        x, y >= 0
* Optimal: x = 3, y = 1, objective 11.
NAME          TINY_LP
OBJSENSE
    MAX
ROWS
 N  profit
 L  c1
 L  c2
COLUMNS
    x         profit       3.0   c1           1.0
    x         c2           1.0
    y         profit       2.0   c1           1.0
    y         c2           3.0
RHS
    rhs       c1           4.0   c2           6.0
BOUNDS
 UP bnd       x            3.0
ENDATA
