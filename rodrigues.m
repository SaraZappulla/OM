function v_rotated = rodrigues(v_initial,u,turning_angle)
    
    v_rotated = v_initial*cos(turning_angle)+cross(u,v_initial)*sin(turning_angle)...
        +dot(u,v_initial).*u.*(1-cos(turning_angle));

end