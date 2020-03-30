function visualize_xnorm(xnorm)
% @brief, plot iteration process
% @param: xnorm, xnorm data

    plot(1:length(xnorm),xnorm);
    xlabel('iterations');
    ylabel('xnorm');
end