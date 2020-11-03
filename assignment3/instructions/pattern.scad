module hexagon(radius = 12, cx = 0, xy = 0, scale = 1)
{
	
	//Points
	p0 = [radius*cos(30), radius*sin(30)];
	p1 = [radius*cos(90), radius*sin(90)];
	p2 = [radius*cos(150), radius*sin(150)];
	p3 = [radius*cos(210), radius*sin(210)];
	p4 = [radius*cos(270), radius*sin(270)];
	p5 = [radius*cos(330), radius*sin(330)];
	
	points = [p0, p1, p2, p3, p4, p5];
	polygon(points);
	}
	
//hexagon();

module planPattern(scalefac = 0.85,R = 360,dx = 36){
	union(){
		dy = dx * 2 * sin(60);
		r = dx/(2 * sin(60))*scalefac;
		//Calculate Grid
		nx = R/dx * 1.2;
		ny = R/dy * 1.2;

		for(i=[(1 - dx*nx)/2:dx:(dx*nx - 1)/2])
			for(j = [(1-dy*ny)/2:dy:(dy*ny - 1)/2])
			{
				//linear_extrude(height = 2*h)
				translate([i, j,0])
					hexagon(radius = r);
				//linear_extrude(height = 2*h)
				translate([i+dx/2, j+dx*sin(60),0])
					hexagon(radius = r);
				}
			}
	}


count = 5;
h = 5;
Rout = 200;
Roffset = 5;
unit = 36;
angle = 20;
	
Rin = Rout - Roffset;

for(i=[0:1:count - 1])
	translate([0,0,h*i])
	linear_extrude(height = h)
	rotate([0,0,angle*i])
	difference(){
		circle(r = Rout);
		intersection()
		{
			//translate([0,0,0.5*h])
			//linear_extrude(height = h)
			circle(r = Rin);
			planPattern(scalefac = 1 - Roffset/unit, R = Rin*2);
			}
		}
		