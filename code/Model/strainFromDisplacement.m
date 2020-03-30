function [e_xx, e_xy, e_yy] = strainFromDisplacement(result, querypoints)    

    [gradx, grady] = evaluateGradient(result,querypoints,[1,2]);
    dudx = gradx(:,1);
    dudy = grady(:,1);
    dvdx = gradx(:,2);
    dvdy = grady(:,2);
    
    e_xx = dudx;
    e_xy = 1/2 * (dudy + dvdx);
    e_yy = dvdy;
end