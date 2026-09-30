NAME          TINY
OBJSENSE
 MIN
ROWS
 N  COST
 G  DEMAND
 L  CAPACITY
COLUMNS
    X1        COST        1          DEMAND      1
    X1        CAPACITY    1
    X2        COST        2          DEMAND      1
    X2        CAPACITY    1
RHS
    RHS1      DEMAND      3          CAPACITY    5
BOUNDS
 UP BND1      X1          4
 UP BND1      X2          4
ENDATA
