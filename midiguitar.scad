//translate([-15, 47.5, pcbO]) import("midiguitarpcb.stl");

w = 68.5;
d = 33;
h = 28;
r = 2;
frontW = 20;
wedgeW = w + r;
wedgeD = 1.5;
wedgeH = 2;
pcbO = 4;

difference() {
	hull() {
		for (x = [0, w]) {
			for (y = [0, d]) {
				for (z = [0, h]) {
					translate([x, y, z]) rotate([0, 45, x > 0 != y > 0 ? -45 : 45]) cube(r, true);
				}
			}
		}
	}
	translate([0, wedgeD, 0]) cube([w + r, d - 2 * wedgeD, h]);
	translate([0, 0, pcbO]) cube([wedgeW, d, wedgeH]);
    translate([w - frontW + r, 0, 0]) cube([frontW, d, h]);
	translate([0, d - 16.5, 7.5 + pcbO + wedgeH]) rotate([0, -90, 0]) cylinder(d = 12, h = r);
}
