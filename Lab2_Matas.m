% Pirmo sluoksnio rysiu svoriai
sl1 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];
%pirmo sluoksnio rezultatai
res1 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];


% Antro sluoksnio rysiu svoriai
sl2 = [randn(1),randn(1),randn(1),randn(1),randn(1),randn(1)];
%antro sluoksnio rezultatas
res2 = randn(1);

n=0.145;

for k = 0:1:100000
    for x = 0.1:1/22:1
        % Apsirašome d
        d = get_d(x);    
        % Pirmo sluoksnio atsakas
        vx_1 = [];
        for i = 1:1:length(sl1)
            vx_1(i)=x*sl1(i)+res1(i);
        end        
        
        % Aktyvioji funkcija
        yx_1 = [];
        for i = 1:1:length(sl1)
            yx_1(i)=1/(1+exp(-vx_1(i)));
        end
    
        % Antro sluoksnio atsakas
        v = 0;
        for i = 1:1:length(sl2)
            v = v + yx_1(i)*sl2(i);
        end

        % klaida
        e=d-v;
    
        % Atnaujinam koficientus
        for i = 1:1:length(sl2)
            sl2(i) = sl2(i)+n*e*yx_1(i);
        end

        res2 = res2 + n * e;
        
        % Pirmo sluoksnio atsakus        
        deltax_1 = [];
        for i = 1:1:length(sl2)
            deltax_1(i)=(yx_1(i)*(1-yx_1(i)))*e*sl2(i);
        end
        
        for i = 1:1:length(sl1)           
            sl1(i)=sl1(i)+n*deltax_1(i)*x;
        end
    
        for i = 1:1:length(deltax_1)           
        res1(i)=res1(i)+n*deltax_1(i);
        end
    end
end

yf = [];
i = 0;

for x = 0.1:1/50:1
    d = get_d(x);

    % Pirmo sluoksio atsakas
    vx_1 = [];
    for i = 1:1:length(sl1)
        vx_1(i)=x*sl1(i)+res1(i);
    end

    % aktyvioji funkcija
    yx_1 = [];
    for i = 1:1:length(vx_1)
        yx_1(i)=1/(1+exp(-vx_1(i)));
    end
    
    % Antro sluoksnio atsakas
    y = 0;
    for i = 1:1:length(sl2)
        y = y + yx_1(i) * sl2(i);
    end
    yf = [yf, y];
end

x = 0.1:1/50:1;
d = get_d(x);

plot ( x, d, '--', x, yf, 'r*')

function d=get_d(x)
    d =  (1 + 0.6 * sin(2*pi*x/0.7)) + 0.3 * sin(2*pi*x) / 2;
end