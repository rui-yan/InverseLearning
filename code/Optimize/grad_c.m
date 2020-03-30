function [gradc_x, gradc_y] = grad_c(x,y, filename, epsilon, stepsize)
    dx = stepsize;
    dy = stepsize;
    c_interp = concentration(filename, epsilon);
    gradc_x = (c_interp(x+dx,y)-c_interp(x-dx,y))./(2*stepsize);
    gradc_y = (c_interp(x,y+dy)-c_interp(x,y-dy))./(2*stepsize);
end

