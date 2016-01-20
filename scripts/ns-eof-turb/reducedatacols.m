function [ col ] = reducedatacols( t )
%% Function returning the meaning of each column (type-sensitive)
%  

    col = struct();

    switch t
        case 'dns'
            col.P      = toD('A');
            col.U1     = toD('B');
            col.U2     = toD('C');
            col.U3     = toD('D');
            col.ARCL   = toD('F');
        case 'aturb'
            col.P      = toD('A');
            col.U1     = toD('B');
            col.U2     = toD('C');
            col.U3     = toD('D');
            
            col.H      = toD('E');
            col.NUT    = toD('F');
            col.LM     = toD('G');
            col.u      = toD('H');
            
            col.ARCL   = toD('J');
        case 'keturb'
            col.P      = toD('A');
            col.U1     = toD('B');
            col.U2     = toD('C');
            col.U3     = toD('D');
            col.H      = toD('E');
            col.NUT    = toD('F');
            col.Pa     = toD('N');
            
            col.TKE    = toD('G');
            col.EPSILON= toD('H');
            col.F1     = toD('I');
            col.F2     = toD('J');
            col.FMU    = toD('K');
            col.SijSij = toD('L');
            col.u      = toD('M');
            
            col.ARCL   = toD('S');
        case 'keturb-de'
            col.P      = toD('A');
            col.U1     = toD('B');
            col.U2     = toD('C');
            col.U3     = toD('D');
            col.H      = toD('E');
            col.NUT    = toD('F');
            col.Pa     = toD('N');
            
            col.TKE    = toD('G');
            col.EPSILON= toD('H');
            col.F1     = toD('I');
            col.F2     = toD('J');
            col.FMU    = toD('K');
            col.SijSij = toD('L');
            col.u      = toD('M');
            
            col.D      = toD('O');
            col.E      = toD('P');
            
            col.ARCL   = toD('U');
    end


    function res = toD(c)
        res = double(c)-double('A')+1;
    end
    
end

