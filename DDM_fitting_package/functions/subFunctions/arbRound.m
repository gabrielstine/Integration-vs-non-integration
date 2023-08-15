function vRounded = arbRound(arbValues, v)

vRounded = interp1(arbValues,arbValues,v,'nearest','extrap') ;
% vRounded = interp1(arbValues,arbValues,v,'previous','extrap') ;

end