function [cone] = cone_list(Ks, threshold, alt)
    %list of cones
    %If the PSD block has size greater than threshold, use the cone in alt
    %otherwise, it's a small block so use psd
    cone = cell(length(Ks), 1);
    for i = 1:length(Ks)
        if Ks(i) <= threshold
            cone{i} = 'psd';
        else
            cone{i} = alt;
        end
    end
end