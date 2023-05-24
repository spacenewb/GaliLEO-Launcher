function plotOrbit(a,e,i,OM,om,th0,thf,dth,mu,color,set,marker)

if th0>thf
    thf=thf+2*pi;
end

if th0==thf
    lwidth=30;
    setwidth='MarkerSize';
else
    lwidth=2;
    setwidth='LineWidth';
end



if set=='3D'
    RR=[];
    for th = th0:dth:thf
        [rr,vv]=par2car(a,e,i,OM,om,th,mu);
        RR=[RR,rr];
    end
    hold on
    grid on

    plot3(RR(1,:),RR(2,:),RR(3,:),color,setwidth,lwidth);

    if marker
        plot3(0,0,0,'.black','MarkerSize',40)
    end
    xlabel('X');
    ylabel('Y');
    zlabel('Z');


elseif set=='2D'
    X=[];
    Y=[];
    i = 0;
    for th = th0:dth:thf
        [ rr_orb, ~] = par2car (a, e, 0, OM, om, th, mu ) ;
        X = [X rr_orb(1)];
        Y = [Y rr_orb(2)];
    end
    hold on
    grid on
    plot(X,Y,color,setwidth,lwidth)
    if marker
        plot(0,0,'.black','MarkerSize',40)
    end
    xlabel('X');
    ylabel('Y');
    
else
    error('tipo plot sbagliato')
end

        
        
    
    





