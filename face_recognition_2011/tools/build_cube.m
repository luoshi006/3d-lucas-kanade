function cube = build_cube(w, h, d)

cube = [
        1 1 1;
        1 h 1;
        w h 1;
        w 1 1;
        1 1 1;
        1 1 d;
        1 h d; 
        1 h 1; 
        1 h d;
        w h d; 
        w h 1; 
        w h d;
        w 1 d; 
        w 1 1; 
        w 1 d;
        1 1 d;
       ]';
end