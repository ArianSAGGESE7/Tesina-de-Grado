function [flip_vec] = flip_ro(vec_ro)
flip_vec = ones(1,length(vec_ro));
    for i=1:length(vec_ro)
        flip_vec(i) = vec_ro(length(vec_ro)-i+1);
    end
end