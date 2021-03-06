function assertMiscEqual( misc, Tmisc )

    assertElementsAlmostEqual(misc.A, Tmisc.A);
    assertElementsAlmostEqual(misc.B, Tmisc.B);
    assertElementsAlmostEqual(misc.intervalSz, Tmisc.intervalSz);
    assertEqual(length(misc.Vm), length(Tmisc.Vm));
    
    for n = 1:length(misc.Vm)
        assertElementsAlmostEqual(misc.Vm{n}, Tmisc.Vm{n});
    end
    
    if ~(all(misc.x(:) + 1 == Tmisc.x'))
        warning('optimal solution did not match; others may fail');
    end    

end

