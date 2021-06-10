classdef mainP < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        secuenciahnEditField            matlab.ui.control.EditField
        secuenciahnEditFieldLabel       matlab.ui.control.Label
        secuenciaxnEditField            matlab.ui.control.EditField
        secuenciaxnEditFieldLabel       matlab.ui.control.Label
        ReflejarButton_2                matlab.ui.control.Button
        FactorInterpolacionEditField_2  matlab.ui.control.NumericEditField
        FactorInterpolacionEditField_2Label  matlab.ui.control.Label
        DiezmacionEditField_2           matlab.ui.control.NumericEditField
        DiezmacionEditField_2Label      matlab.ui.control.Label
        AtenuacionAmplificacionEditField_2  matlab.ui.control.NumericEditField
        AtenuacionAmplificacionEditField_2Label  matlab.ui.control.Label
        DesplazamientoEditField_2       matlab.ui.control.NumericEditField
        DesplazamientoEditField_2Label  matlab.ui.control.Label
        LinealButton_2                  matlab.ui.control.Button
        EscalonButton_3                 matlab.ui.control.Button
        CerosButton_2                   matlab.ui.control.Button
        DiezmarButton_2                 matlab.ui.control.Button
        AtenuarAmplificarButton_2       matlab.ui.control.Button
        DesplazarButton_2               matlab.ui.control.Button
        OrigenhnEditField               matlab.ui.control.NumericEditField
        OrigenhnLabel                   matlab.ui.control.Label
        OrigenxnEditField               matlab.ui.control.NumericEditField
        OrigenxnEditFieldLabel          matlab.ui.control.Label
        LinealButton                    matlab.ui.control.Button
        EscalonButton_2                 matlab.ui.control.Button
        SecuenciasButton                matlab.ui.control.Button
        ConvolucionarButton             matlab.ui.control.Button
        RestarButton                    matlab.ui.control.Button
        SumarButton                     matlab.ui.control.Button
        CerosButton                     matlab.ui.control.Button
        DiezmarButton                   matlab.ui.control.Button
        ReflejarButton                  matlab.ui.control.Button
        AtenuarAmplificarButton         matlab.ui.control.Button
        DesplazarButton                 matlab.ui.control.Button
        FactorInterpolacionEditField    matlab.ui.control.NumericEditField
        FactorInterpolacionEditFieldLabel  matlab.ui.control.Label
        DiezmacionEditField             matlab.ui.control.NumericEditField
        DiezmacionEditFieldLabel        matlab.ui.control.Label
        AtenuacionAmplificacionEditField  matlab.ui.control.NumericEditField
        AtenuacionAmplificacionEditFieldLabel  matlab.ui.control.Label
        DesplazamientoEditField         matlab.ui.control.NumericEditField
        DesplazamientoEditFieldLabel    matlab.ui.control.Label
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes_4                        matlab.ui.control.UIAxes
        UIAxes_3                        matlab.ui.control.UIAxes
        UIAxes_2                        matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        X
        LongitudX
        OrigenX
        InferiorX
        SuperiorX
        Y
        X2
        LongitudX2
        OrigenX2
        InferiorX2
        SuperiorX2
        Y2
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: SecuenciasButton
        function cargaSecuencia(app, event)
% load Data.dat;
% app.X = Data(1,:);
app.X=[];
a=app.secuenciaxnEditField.Value;
ac=[];
b=erase(a,{'[',']'});
a=split(b,{' ',','});
ac=str2double(a);
for z=1:length(ac)
    app.X(z)=ac(z);
end
app.LongitudX = length(app.X);
app.OrigenX = app.OrigenxnEditField.Value ;
app.InferiorX = (-1) * ( app.OrigenX - 1) ;
app.SuperiorX = app.LongitudX - app.OrigenX ;
app.Y = app.InferiorX : app.SuperiorX ;
            stem(app.UIAxes,app.Y,app.X);

% load Data2.dat;
% app.X2 = Data2(1,:);
app.X2=[];
a2=app.secuenciahnEditField.Value;
ac2=[];
b2=erase(a2,{'[',']'});
a2=split(b2,{' ',','});
ac2=str2double(a2);
for z=1:length(ac2)
    app.X2(z)=ac2(z);
end
app.LongitudX2 = length(app.X2);
app.OrigenX2 = app.OrigenhnEditField.Value ;
app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
            stem(app.UIAxes_3,app.Y2,app.X2);
        end

        % Button pushed function: DesplazarButton
        function DesplazarButtonPushed(app, event)
Desplazamiento = app.DesplazamientoEditField.Value ;
app.OrigenX = app.OrigenX + Desplazamiento ;
app.InferiorX = (-1) * ( app.OrigenX - 1) ;
app.SuperiorX = app.LongitudX - app.OrigenX ;
if app.InferiorX > 0
    E = [ 0 ] ;
    E = repmat(E,1, app.InferiorX );    
    app.X =  cat( 2 , E , app.X  );    
    app.LongitudX = length( app.X );    
    app.OrigenX = app.OrigenX + app.InferiorX ;
    app.InferiorX = (-1) * ( app.OrigenX - 1) ;
    app.SuperiorX = app.LongitudX - app.OrigenX ;
end
if app.SuperiorX < 0
    E = [ 0 ] ;
    E = repmat(E,1, app.SuperiorX*(-1) );    
    app.X =  cat( 2 , app.X , E );    
    app.LongitudX = length( app.X );
    app.SuperiorX = app.LongitudX - app.OrigenX ;
end
app.Y = app.InferiorX : app.SuperiorX ;
stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: AtenuarAmplificarButton
        function AtenuarAmplificarButtonPushed(app, event)
app.X = app.AtenuacionAmplificacionEditField.Value * app.X;
app.InferiorX = (-1) * ( app.OrigenX - 1) ;
app.SuperiorX = app.LongitudX - app.OrigenX ;
app.Y = app.InferiorX : app.SuperiorX ;
stem(app.UIAxes_2,app.Y,app.X);            
        end

        % Button pushed function: ReflejarButton
        function ReflejarButtonPushed(app, event)
    Aux = app.InferiorX ;
    app.InferiorX = app.SuperiorX * ( -1 ) ;
    app.SuperiorX = Aux * ( -1 ) ;
    app.X = fliplr(app.X) ;
    app.LongitudX = length(app.X);
    app.OrigenX = ( app.LongitudX - app.OrigenX ) + 1 ;
    app.Y = app.InferiorX : app.SuperiorX ;
    stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: DiezmarButton
        function DiezmarButtonPushed(app, event)
            factor = app.DiezmacionEditField.Value ;

Origen = app.OrigenX ;
X1 = app.X(Origen:factor:end) ;
X1 = X1( 2 : end ) ;
app.X = fliplr(app.X) ;
Longitud = length(app.X);
Origen = ( Longitud - Origen ) + 1 ;
Xa2 = app.X(Origen:factor:end) ;
Xa2 = Xa2( 2 : end ) ;
Xa2 = fliplr(Xa2) ;
R = cat( 2, Xa2, app.X(Origen) );
R = cat( 2, R, X1 );
app.X = R ;
app.LongitudX = length(R);
app.OrigenX = length(Xa2) + 1 ;
app.InferiorX = (-1) * ( app.OrigenX - 1) ;
app.SuperiorX = app.LongitudX - app.OrigenX ;
app.Y = app.InferiorX : app.SuperiorX ;
stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: CerosButton
        function CerosButtonPushed(app, event)
        E = [ 0 ] ;
        E = repmat(E,1, app.FactorInterpolacionEditField.Value );        
        Origen = app.OrigenX ;        
        R1 = [] ;
        for i=Origen : app.LongitudX %Ciclo con el cual se llenará el arreglo
            if i==Origen
                
            else
                Aux = cat( 2 , E , app.X(i) );
                R1 = cat( 2 , R1 , Aux );
            end
        end        
        app.X = fliplr(app.X) ;
        Longitud = length(app.X);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=Origen : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
                
            else
                Aux = cat( 2 , E , app.X(i) );
                R2 = cat( 2 , R2 , Aux );
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X(Origen) );
        R = cat( 2, R, R1 );
        app.X = R ;
        app.LongitudX = length(R);
        app.OrigenX = length(R2) + 1 ;
        app.InferiorX = (-1) * ( app.OrigenX - 1) ;
        app.SuperiorX = app.LongitudX - app.OrigenX ;
        app.Y = app.InferiorX : app.SuperiorX ;
        stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: EscalonButton_2
        function EscalonButton_2Pushed(app, event)
        Origen = app.OrigenX ;        
        R1 = [] ;
        for i=Origen : app.LongitudX %Ciclo con el cual se llenará el arreglo
            if i==Origen
                E = [ app.X(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField.Value );
            else
                E = [ app.X(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField.Value );
                Aux = cat( 2 , E , app.X(i) );
                R1 = cat( 2 , R1 , Aux );
            end
        end        
        app.X = fliplr(app.X) ;
        Longitud = length(app.X);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=Origen : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
                E = [ app.X(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField.Value );
            else
                E = [ app.X(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField.Value );
                Aux = cat( 2 , E , app.X(i) );
                R2 = cat( 2 , R2 , Aux );
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X(Origen) );
        R = cat( 2, R, R1 );
        app.X = R ;
        app.LongitudX = length(R);
        app.OrigenX = length(R2) + 1 ;
        app.InferiorX = (-1) * ( app.OrigenX - 1) ;
        app.SuperiorX = app.LongitudX - app.OrigenX ;       
        app.Y = app.InferiorX : app.SuperiorX ;
        stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: LinealButton
        function LinealButtonPushed(app, event)
        Origen = app.OrigenX ;        
        R1 = [] ;
        for i=Origen : app.LongitudX %Ciclo con el cual se llenará el arreglo
            if i==Origen
            else
                if app.X(i-1) < app.X(i)
                    E = linspace( app.X(i-1) , app.X(i) , app.FactorInterpolacionEditField.Value+2 ) ;
                    E = E( 2 : end );
                    R1 = cat( 2 , R1 , E );
                else
                    E = linspace( app.X(i) , app.X(i-1) , app.FactorInterpolacionEditField.Value+2 ) ;
                    E = fliplr(E) ;
                    E = E( 2 : end );
                    R1 = cat( 2 , R1 , E );
                end
            end
        end
        app.X = fliplr(app.X) ;
        Longitud = length(app.X);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=Origen : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
            else
                if app.X(i-1) < app.X(i)
                    E = linspace( app.X(i-1) , app.X(i) , app.FactorInterpolacionEditField.Value+2 ) ;
                    E = E( 2 : end );
                    R2 = cat( 2 , R2 , E );
                else
                    E = linspace( app.X(i) , app.X(i-1) , app.FactorInterpolacionEditField.Value+2 ) ;
                    E = fliplr(E) ;
                    E = E( 2 : end );
                    R2 = cat( 2 , R2 , E );
                end
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X(Origen) );
        R = cat( 2, R, R1 );
        app.X = R ;
        app.LongitudX = length(R);
        app.OrigenX = length(R2) + 1 ;
        app.InferiorX = (-1) * ( app.OrigenX - 1) ;
        app.SuperiorX = app.LongitudX - app.OrigenX ;        
        app.Y = app.InferiorX : app.SuperiorX ;
        stem(app.UIAxes_2,app.Y,app.X);
        end

        % Button pushed function: ConvolucionarButton
        function ConvolucionarButtonPushed(app, event)
e = app.X;
h = app.X2;
m = length(h)+length(e)-1;
n = length(h);
w = zeros(m,n);
y = zeros(1,n);
a = 1;
for k=1:m
    for b=n:-1:2
        y(b) = y(b-1);
    end
    if k<=(m-n+1)
        y(1)=e(k);
    else
        y(1) = 0;
    end
        for i=1:n
            w(k,i) = h(i)*y(i);
            a=a+1;
        end
            a=1;
end
    for j=1:m
        xnuevo=0;
        xant=0;
        for c=1:n
            xant=w(j,c);
            xnuevo = xant+xnuevo;
        end
         W(j)=xnuevo;
    end
Longitud= length(W);
Origen =1 + abs(app.InferiorX-app.InferiorX2);
Inferior = (-1) * ( Origen - 1) ;
Superior = Longitud - Origen ;
Ya = Inferior : Superior ;
stem(app.UIAxes2,Ya,W)
        end

        % Button pushed function: SumarButton
        function SumarButtonPushed(app, event)
            arr1=app.X;
            arr2=app.X2;
            inf1=app.InferiorX;
            inf2=app.InferiorX2;
            i=1;
            j=1;
            if inf1 < inf2
                for a=1:length(arr1)
                    R(a)=arr1(a);
                    if inf1<inf2
                        i=i+1;
                    end
                    inf1=inf1+1;
                end
                for b=i:length(arr1)
                    R(b)=arr1(b)+arr2(j);
                    j=j+1;
                    i=i+1;
                end
                for c=j:length(arr2)
                    R(i)=arr2(c);
                    i=i+1;
                end
                Longitud= length(R);
Origen =app.OrigenX;
Inferior = (-1) * ( Origen - 1) ;
Superior = Longitud - Origen ;
Ya = Inferior : Superior ;
stem(app.UIAxes2,Ya,R)
            elseif inf1 >inf2
                for a=1:length(arr2)
                    R(a)=arr2(a);
                    if inf1>inf2
                        i=i+1;
                    end
                    inf2=inf2+1;
                end
                for b=i:length(arr2)
                    R(b)=arr1(j)+arr2(b);
                    j=j+1;
                    i=i+1;
                end
                for c=j:length(arr1)
                    R(i)=arr1(c);
                    i=i+1;
                end
Longitud= length(R);
Origen =app.OrigenX2;
Inferior = (-1) * ( Origen - 1) ;
Superior = Longitud - Origen ;
Ya = Inferior : Superior ;
stem(app.UIAxes2,Ya,R)
            end
        end

        % Button pushed function: RestarButton
        function RestarButtonPushed(app, event)
            arr1=app.X;
            arr2=app.X2;
            inf1=app.InferiorX;
            inf2=app.InferiorX2;
            i=1;
            j=1;
            if inf1 < inf2
                for a=1:length(arr1)
                    R(a)=arr1(a);
                    if inf1<inf2
                        i=i+1;
                    end
                    inf1=inf1+1;
                end
                for b=i:length(arr1)
                    R(b)=arr1(b)-arr2(j);
                    j=j+1;
                    i=i+1;
                end
                for c=j:length(arr2)
                    R(i)=-arr2(c);
                    i=i+1;
                end
                Longitud= length(R);
Origen =app.OrigenX;
Inferior = (-1) * ( Origen - 1) ;
Superior = Longitud - Origen ;
Ya = Inferior : Superior ;
stem(app.UIAxes2,Ya,R)
            elseif inf1 >inf2
                for a=1:length(arr2)
                    R(a)=-arr2(a);
                    if inf1>inf2
                        i=i+1;
                    end
                    inf2=inf2+1;
                end
                for b=i:length(arr2)
                    R(b)=arr1(j)-arr2(b);
                    j=j+1;
                    i=i+1;
                end
                for c=j:length(arr1)
                    R(i)=arr1(c);
                    i=i+1;
                end
Longitud= length(R);
Origen =app.OrigenX2;
Inferior = (-1) * ( Origen - 1) ;
Superior = Longitud - Origen ;
Ya = Inferior : Superior ;
stem(app.UIAxes2,Ya,R)
            end
        end

        % Button pushed function: ReflejarButton_2
        function ReflejarButton_2Pushed(app, event)
            Aux2 = app.InferiorX2 ;
    app.InferiorX2= app.SuperiorX2 * ( -1 ) ;
    app.SuperiorX2 = Aux2 * ( -1 ) ;
    app.X2 = fliplr(app.X2) ;
    app.LongitudX2 = length(app.X2);
    app.OrigenX2 = ( app.LongitudX2 - app.OrigenX2 ) + 1 ;
    app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
    stem(app.UIAxes_4,app.Y2,app.X2);        
        end

        % Button pushed function: AtenuarAmplificarButton_2
        function AtenuarAmplificarButton_2Pushed(app, event)
app.X2 = app.AtenuacionAmplificacionEditField_2.Value * app.X2;
app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
stem(app.UIAxes_4,app.Y2,app.X2);
        end

        % Button pushed function: DesplazarButton_2
        function DesplazarButton_2Pushed(app, event)
            Desplazamiento = app.DesplazamientoEditField_2.Value ;
app.OrigenX2 = app.OrigenX2 + Desplazamiento ;
app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
if app.InferiorX2 > 0
    E = [ 0 ] ;
    E = repmat(E,1, app.InferiorX2 );    
    app.X2 =  cat( 2 , E , app.X2  );    
    app.LongitudX2 = length( app.X2 );    
    app.OrigenX2 = app.OrigenX2 + app.InferiorX2 ;
    app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
    app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
end
if app.SuperiorX2 < 0
    E = [ 0 ] ;
    E = repmat(E,1, app.SuperiorX2*(-1) );    
    app.X2 =  cat( 2 , app.X2 , E );    
    app.LongitudX2 = length( app.X2 );
    app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
end
app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
stem(app.UIAxes_4,app.Y2,app.X2);
        end

        % Button pushed function: DiezmarButton_2
        function DiezmarButton_2Pushed(app, event)
            Origen22 = app.OrigenX2 ;
            factor = app.DiezmacionEditField_2.Value ;
X12 = app.X2(Origen22:factor:end) ;
X12 = X12( 2 : end ) ;
app.X2 = fliplr(app.X2) ;
Longitud22 = length(app.X2);
Origen22 = ( Longitud22 - Origen22 ) + 1 ;
Xa22 = app.X2(Origen22:factor:end) ;
Xa22 = Xa22( 2 : end ) ;
Xa22 = fliplr(Xa22) ;
R2 = cat( 2, Xa22, app.X2(Origen22) );
R2 = cat( 2, R2, X12 );
app.X2 = R2 ;
app.LongitudX2 = length(R2);
app.OrigenX2 = length(Xa22) + 1 ;
app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
stem(app.UIAxes_4,app.Y2,app.X2);
        end

        % Button pushed function: CerosButton_2
        function CerosButton_2Pushed(app, event)
            E = [ 0 ] ;
        E = repmat(E,1, app.FactorInterpolacionEditField_2.Value );        
        Origen = app.OrigenX2 ;        
        R1 = [] ;
        for i=Origen : app.LongitudX2 %Ciclo con el cual se llenará el arreglo
            if i==Origen
                
            else
                Aux = cat( 2 , E , app.X2(i) );
                R1 = cat( 2 , R1 , Aux );
            end
        end        
        app.X2 = fliplr(app.X2) ;
        Longitud = length(app.X2);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=Origen : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
                
            else
                Aux = cat( 2 , E , app.X2(i) );
                R2 = cat( 2 , R2 , Aux );
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X2(Origen) );
        R = cat( 2, R, R1 );
        app.X2 = R ;
        app.LongitudX2 = length(R);
        app.OrigenX2 = length(R2) + 1 ;
        app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
        app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;
        app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
        stem(app.UIAxes_4,app.Y2,app.X2);
        end

        % Button pushed function: EscalonButton_3
        function EscalonButton_3Pushed(app, event)
            Origen = app.OrigenX2 ;        
        R1 = [] ;
        for i=Origen : app.LongitudX2 %Ciclo con el cual se llenará el arreglo
            if i==Origen
                E = [ app.X2(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField_2.Value );
            else
                E = [ app.X2(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField_2.Value );
                Aux = cat( 2 , E , app.X2(i) );
                R1 = cat( 2 , R1 , Aux );
            end
        end        
        app.X2 = fliplr(app.X2) ;
        Longitud = length(app.X2);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=Origen : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
                E = [ app.X2(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField_2.Value );
            else
                E = [ app.X2(i) ] ;
                E = repmat(E,1, app.FactorInterpolacionEditField_2.Value );
                Aux = cat( 2 , E , app.X2(i) );
                R2 = cat( 2 , R2 , Aux );
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X2(Origen) );
        R = cat( 2, R, R1 );
        app.X2 = R ;
        app.LongitudX2 = length(R);
        app.OrigenX2 = length(R2) + 1 ;
        app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
        app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;       
        app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
        stem(app.UIAxes_4,app.Y2,app.X2);
        end

        % Button pushed function: LinealButton_2
        function LinealButton_2Pushed(app, event)
            Origen = app.OrigenX2 ;        
        R1 = [] ;
        for i=Origen : app.LongitudX2 %Ciclo con el cual se llenará el arreglo
            if i==Origen
            else
                if app.X2(i-1) < app.X2(i)
                    E = linspace( app.X2(i-1) , app.X2(i) , app.FactorInterpolacionEditField_2.Value+2 ) ;
                    E = E( 2 : end );
                    R1 = cat( 2 , R1 , E );
                else
                    E = linspace( app.X2(i) , app.X2(i-1) , app.FactorInterpolacionEditField_2.Value+2 ) ;
                    E = fliplr(E) ;
                    E = E( 2 : end );
                    R1 = cat( 2 , R1 , E );
                end
            end
        end
        app.X2 = fliplr(app.X2) ;
        Longitud = length(app.X2);
        Origen = ( Longitud - Origen ) + 1 ;        
        R2 = [] ;
        for i=abs(Origen) : Longitud %Ciclo con el cual se llenará el arreglo
            if i==Origen
            else
                if app.X2(i-1) < app.X2(i)
                    E = linspace( app.X2(i-1) , app.X2(i) , app.FactorInterpolacionEditField_2.Value+2 ) ;
                    E = E( 2 : end );
                    R2 = cat( 2 , R2 , E );
                else
                    E = linspace( app.X2(i) , app.X2(i-1) , app.FactorInterpolacionEditField_2.Value+2 ) ;
                    E = fliplr(E) ;
                    E = E( 2 : end );
                    R2 = cat( 2 , R2 , E );
                end
            end
        end        
        R2 = fliplr(R2) ;        
        R = cat( 2, R2, app.X2(Origen) );
        R = cat( 2, R, R1 );
        app.X2 = R ;
        app.LongitudX2 = length(R);
        app.OrigenX2 = length(R2) + 1 ;
        app.InferiorX2 = (-1) * ( app.OrigenX2 - 1) ;
        app.SuperiorX2 = app.LongitudX2 - app.OrigenX2 ;        
        app.Y2 = app.InferiorX2 : app.SuperiorX2 ;
        stem(app.UIAxes_4,app.Y2,app.X2);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1264 734];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'x(n)')
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [16 487 434 234];

            % Create UIAxes_2
            app.UIAxes_2 = uiaxes(app.UIFigure);
            title(app.UIAxes_2, 'Resultado x(n)')
            xlabel(app.UIAxes_2, 'X')
            ylabel(app.UIAxes_2, 'Y')
            zlabel(app.UIAxes_2, 'Z')
            app.UIAxes_2.Position = [16 257 434 231];

            % Create UIAxes_3
            app.UIAxes_3 = uiaxes(app.UIFigure);
            title(app.UIAxes_3, 'h(n)')
            xlabel(app.UIAxes_3, 'X')
            ylabel(app.UIAxes_3, 'Y')
            zlabel(app.UIAxes_3, 'Z')
            app.UIAxes_3.Position = [834 487 418 234];

            % Create UIAxes_4
            app.UIAxes_4 = uiaxes(app.UIFigure);
            title(app.UIAxes_4, 'Resultado h(n)')
            xlabel(app.UIAxes_4, 'X')
            ylabel(app.UIAxes_4, 'Y')
            zlabel(app.UIAxes_4, 'Z')
            app.UIAxes_4.Position = [839 257 413 231];

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.UIFigure);
            title(app.UIAxes2, 'x(n) y h(n)')
            xlabel(app.UIAxes2, 'X')
            ylabel(app.UIAxes2, 'Y')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.Position = [400 21 466 228];

            % Create DesplazamientoEditFieldLabel
            app.DesplazamientoEditFieldLabel = uilabel(app.UIFigure);
            app.DesplazamientoEditFieldLabel.HorizontalAlignment = 'right';
            app.DesplazamientoEditFieldLabel.Position = [61 187 92 22];
            app.DesplazamientoEditFieldLabel.Text = 'Desplazamiento';

            % Create DesplazamientoEditField
            app.DesplazamientoEditField = uieditfield(app.UIFigure, 'numeric');
            app.DesplazamientoEditField.Position = [168 186 32 24];
            app.DesplazamientoEditField.Value = 1;

            % Create AtenuacionAmplificacionEditFieldLabel
            app.AtenuacionAmplificacionEditFieldLabel = uilabel(app.UIFigure);
            app.AtenuacionAmplificacionEditFieldLabel.HorizontalAlignment = 'right';
            app.AtenuacionAmplificacionEditFieldLabel.Position = [4 141 140 22];
            app.AtenuacionAmplificacionEditFieldLabel.Text = 'Atenuacion/Amplificacion';

            % Create AtenuacionAmplificacionEditField
            app.AtenuacionAmplificacionEditField = uieditfield(app.UIFigure, 'numeric');
            app.AtenuacionAmplificacionEditField.Position = [167 138 32 24];
            app.AtenuacionAmplificacionEditField.Value = 2;

            % Create DiezmacionEditFieldLabel
            app.DiezmacionEditFieldLabel = uilabel(app.UIFigure);
            app.DiezmacionEditFieldLabel.HorizontalAlignment = 'right';
            app.DiezmacionEditFieldLabel.Position = [84 95 68 22];
            app.DiezmacionEditFieldLabel.Text = 'Diezmacion';

            % Create DiezmacionEditField
            app.DiezmacionEditField = uieditfield(app.UIFigure, 'numeric');
            app.DiezmacionEditField.Position = [167 94 32 24];
            app.DiezmacionEditField.Value = 2;

            % Create FactorInterpolacionEditFieldLabel
            app.FactorInterpolacionEditFieldLabel = uilabel(app.UIFigure);
            app.FactorInterpolacionEditFieldLabel.HorizontalAlignment = 'right';
            app.FactorInterpolacionEditFieldLabel.Position = [41 57 112 22];
            app.FactorInterpolacionEditFieldLabel.Text = 'Factor Interpolacion';

            % Create FactorInterpolacionEditField
            app.FactorInterpolacionEditField = uieditfield(app.UIFigure, 'numeric');
            app.FactorInterpolacionEditField.Position = [168 56 32 24];
            app.FactorInterpolacionEditField.Value = 2;

            % Create DesplazarButton
            app.DesplazarButton = uibutton(app.UIFigure, 'push');
            app.DesplazarButton.ButtonPushedFcn = createCallbackFcn(app, @DesplazarButtonPushed, true);
            app.DesplazarButton.Position = [219 187 100 22];
            app.DesplazarButton.Text = 'Desplazar';

            % Create AtenuarAmplificarButton
            app.AtenuarAmplificarButton = uibutton(app.UIFigure, 'push');
            app.AtenuarAmplificarButton.ButtonPushedFcn = createCallbackFcn(app, @AtenuarAmplificarButtonPushed, true);
            app.AtenuarAmplificarButton.Position = [205 139 114 22];
            app.AtenuarAmplificarButton.Text = 'Atenuar/Amplificar';

            % Create ReflejarButton
            app.ReflejarButton = uibutton(app.UIFigure, 'push');
            app.ReflejarButton.ButtonPushedFcn = createCallbackFcn(app, @ReflejarButtonPushed, true);
            app.ReflejarButton.Position = [133 227 100 22];
            app.ReflejarButton.Text = 'Reflejar';

            % Create DiezmarButton
            app.DiezmarButton = uibutton(app.UIFigure, 'push');
            app.DiezmarButton.ButtonPushedFcn = createCallbackFcn(app, @DiezmarButtonPushed, true);
            app.DiezmarButton.Position = [219 95 100 22];
            app.DiezmarButton.Text = 'Diezmar';

            % Create CerosButton
            app.CerosButton = uibutton(app.UIFigure, 'push');
            app.CerosButton.ButtonPushedFcn = createCallbackFcn(app, @CerosButtonPushed, true);
            app.CerosButton.Position = [219 57 100 22];
            app.CerosButton.Text = 'Ceros';

            % Create SumarButton
            app.SumarButton = uibutton(app.UIFigure, 'push');
            app.SumarButton.ButtonPushedFcn = createCallbackFcn(app, @SumarButtonPushed, true);
            app.SumarButton.Position = [460 307 100 22];
            app.SumarButton.Text = 'Sumar';

            % Create RestarButton
            app.RestarButton = uibutton(app.UIFigure, 'push');
            app.RestarButton.ButtonPushedFcn = createCallbackFcn(app, @RestarButtonPushed, true);
            app.RestarButton.Position = [593 307 100 22];
            app.RestarButton.Text = 'Restar';

            % Create ConvolucionarButton
            app.ConvolucionarButton = uibutton(app.UIFigure, 'push');
            app.ConvolucionarButton.ButtonPushedFcn = createCallbackFcn(app, @ConvolucionarButtonPushed, true);
            app.ConvolucionarButton.Position = [726 307 100 22];
            app.ConvolucionarButton.Text = 'Convolucionar';

            % Create SecuenciasButton
            app.SecuenciasButton = uibutton(app.UIFigure, 'push');
            app.SecuenciasButton.ButtonPushedFcn = createCallbackFcn(app, @cargaSecuencia, true);
            app.SecuenciasButton.Position = [583 487 100 22];
            app.SecuenciasButton.Text = 'Secuencias';

            % Create EscalonButton_2
            app.EscalonButton_2 = uibutton(app.UIFigure, 'push');
            app.EscalonButton_2.ButtonPushedFcn = createCallbackFcn(app, @EscalonButton_2Pushed, true);
            app.EscalonButton_2.Position = [71 20 100 22];
            app.EscalonButton_2.Text = 'Escalon';

            % Create LinealButton
            app.LinealButton = uibutton(app.UIFigure, 'push');
            app.LinealButton.ButtonPushedFcn = createCallbackFcn(app, @LinealButtonPushed, true);
            app.LinealButton.Position = [219 20 100 22];
            app.LinealButton.Text = 'Lineal';

            % Create OrigenxnEditFieldLabel
            app.OrigenxnEditFieldLabel = uilabel(app.UIFigure);
            app.OrigenxnEditFieldLabel.HorizontalAlignment = 'right';
            app.OrigenxnEditFieldLabel.Position = [460 622 66 22];
            app.OrigenxnEditFieldLabel.Text = 'Origen x(n)';

            % Create OrigenxnEditField
            app.OrigenxnEditField = uieditfield(app.UIFigure, 'numeric');
            app.OrigenxnEditField.Position = [573 622 100 22];
            app.OrigenxnEditField.Value = 3;

            % Create OrigenhnLabel
            app.OrigenhnLabel = uilabel(app.UIFigure);
            app.OrigenhnLabel.HorizontalAlignment = 'right';
            app.OrigenhnLabel.Position = [760 545 66 22];
            app.OrigenhnLabel.Text = 'Origen h(n)';

            % Create OrigenhnEditField
            app.OrigenhnEditField = uieditfield(app.UIFigure, 'numeric');
            app.OrigenhnEditField.Position = [613 545 100 22];
            app.OrigenhnEditField.Value = 2;

            % Create DesplazarButton_2
            app.DesplazarButton_2 = uibutton(app.UIFigure, 'push');
            app.DesplazarButton_2.ButtonPushedFcn = createCallbackFcn(app, @DesplazarButton_2Pushed, true);
            app.DesplazarButton_2.Position = [1094 181 100 22];
            app.DesplazarButton_2.Text = 'Desplazar';

            % Create AtenuarAmplificarButton_2
            app.AtenuarAmplificarButton_2 = uibutton(app.UIFigure, 'push');
            app.AtenuarAmplificarButton_2.ButtonPushedFcn = createCallbackFcn(app, @AtenuarAmplificarButton_2Pushed, true);
            app.AtenuarAmplificarButton_2.Position = [1080 133 114 22];
            app.AtenuarAmplificarButton_2.Text = 'Atenuar/Amplificar';

            % Create DiezmarButton_2
            app.DiezmarButton_2 = uibutton(app.UIFigure, 'push');
            app.DiezmarButton_2.ButtonPushedFcn = createCallbackFcn(app, @DiezmarButton_2Pushed, true);
            app.DiezmarButton_2.Position = [1094 89 100 22];
            app.DiezmarButton_2.Text = 'Diezmar';

            % Create CerosButton_2
            app.CerosButton_2 = uibutton(app.UIFigure, 'push');
            app.CerosButton_2.ButtonPushedFcn = createCallbackFcn(app, @CerosButton_2Pushed, true);
            app.CerosButton_2.Position = [1094 51 100 22];
            app.CerosButton_2.Text = 'Ceros';

            % Create EscalonButton_3
            app.EscalonButton_3 = uibutton(app.UIFigure, 'push');
            app.EscalonButton_3.ButtonPushedFcn = createCallbackFcn(app, @EscalonButton_3Pushed, true);
            app.EscalonButton_3.Position = [946 14 100 22];
            app.EscalonButton_3.Text = 'Escalon';

            % Create LinealButton_2
            app.LinealButton_2 = uibutton(app.UIFigure, 'push');
            app.LinealButton_2.ButtonPushedFcn = createCallbackFcn(app, @LinealButton_2Pushed, true);
            app.LinealButton_2.Position = [1094 14 100 22];
            app.LinealButton_2.Text = 'Lineal';

            % Create DesplazamientoEditField_2Label
            app.DesplazamientoEditField_2Label = uilabel(app.UIFigure);
            app.DesplazamientoEditField_2Label.HorizontalAlignment = 'right';
            app.DesplazamientoEditField_2Label.Position = [936 181 92 22];
            app.DesplazamientoEditField_2Label.Text = 'Desplazamiento';

            % Create DesplazamientoEditField_2
            app.DesplazamientoEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.DesplazamientoEditField_2.Position = [1043 180 32 24];
            app.DesplazamientoEditField_2.Value = 1;

            % Create AtenuacionAmplificacionEditField_2Label
            app.AtenuacionAmplificacionEditField_2Label = uilabel(app.UIFigure);
            app.AtenuacionAmplificacionEditField_2Label.HorizontalAlignment = 'right';
            app.AtenuacionAmplificacionEditField_2Label.Position = [887 146 140 22];
            app.AtenuacionAmplificacionEditField_2Label.Text = 'Atenuacion/Amplificacion';

            % Create AtenuacionAmplificacionEditField_2
            app.AtenuacionAmplificacionEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.AtenuacionAmplificacionEditField_2.Position = [1042 132 32 24];
            app.AtenuacionAmplificacionEditField_2.Value = 2;

            % Create DiezmacionEditField_2Label
            app.DiezmacionEditField_2Label = uilabel(app.UIFigure);
            app.DiezmacionEditField_2Label.HorizontalAlignment = 'right';
            app.DiezmacionEditField_2Label.Position = [959 89 68 22];
            app.DiezmacionEditField_2Label.Text = 'Diezmacion';

            % Create DiezmacionEditField_2
            app.DiezmacionEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.DiezmacionEditField_2.Position = [1042 88 32 24];
            app.DiezmacionEditField_2.Value = 2;

            % Create FactorInterpolacionEditField_2Label
            app.FactorInterpolacionEditField_2Label = uilabel(app.UIFigure);
            app.FactorInterpolacionEditField_2Label.HorizontalAlignment = 'right';
            app.FactorInterpolacionEditField_2Label.Position = [916 51 112 22];
            app.FactorInterpolacionEditField_2Label.Text = 'Factor Interpolacion';

            % Create FactorInterpolacionEditField_2
            app.FactorInterpolacionEditField_2 = uieditfield(app.UIFigure, 'numeric');
            app.FactorInterpolacionEditField_2.Position = [1043 50 32 24];
            app.FactorInterpolacionEditField_2.Value = 2;

            % Create ReflejarButton_2
            app.ReflejarButton_2 = uibutton(app.UIFigure, 'push');
            app.ReflejarButton_2.ButtonPushedFcn = createCallbackFcn(app, @ReflejarButton_2Pushed, true);
            app.ReflejarButton_2.Position = [1009 219 100 22];
            app.ReflejarButton_2.Text = 'Reflejar';

            % Create secuenciaxnEditFieldLabel
            app.secuenciaxnEditFieldLabel = uilabel(app.UIFigure);
            app.secuenciaxnEditFieldLabel.HorizontalAlignment = 'right';
            app.secuenciaxnEditFieldLabel.Position = [460 663 84 22];
            app.secuenciaxnEditFieldLabel.Text = 'secuencia x(n)';

            % Create secuenciaxnEditField
            app.secuenciaxnEditField = uieditfield(app.UIFigure, 'text');
            app.secuenciaxnEditField.Position = [573 663 253 22];
            app.secuenciaxnEditField.Value = '[1,2,3,4]';

            % Create secuenciahnEditFieldLabel
            app.secuenciahnEditFieldLabel = uilabel(app.UIFigure);
            app.secuenciahnEditFieldLabel.HorizontalAlignment = 'right';
            app.secuenciahnEditFieldLabel.Position = [742 581 84 22];
            app.secuenciahnEditFieldLabel.Text = 'secuencia h(n)';

            % Create secuenciahnEditField
            app.secuenciahnEditField = uieditfield(app.UIFigure, 'text');
            app.secuenciahnEditField.Position = [460 581 253 22];
            app.secuenciahnEditField.Value = '[4,5,6,7]';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = mainP

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end