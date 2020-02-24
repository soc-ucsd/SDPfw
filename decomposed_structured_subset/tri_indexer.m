function [Mi, L] = tri_indexer(num_var)
    %create a moment matrix with indices in Mi and L variables
    m_size = num_var;
    L = m_size*(m_size+1)/2;

    Mi = zeros(m_size);
    %t = 1:m_size;

    %k_first = 1;
    %k_last = num_var;
    k = 0;
    for i = 1:m_size
        k_range = k + (1:(m_size + 1 - i));
        Mi(i, i:m_size) = k_range;
        Mi(i:m_size, i) = k_range;
        k = k + (m_size + 1 - i);
    end


end