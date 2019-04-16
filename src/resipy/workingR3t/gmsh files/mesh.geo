cl1=30;
cl2=5;
cl3=1;


// electrodes in centre line
Point(101) = {0, 0, 0, cl3};
Point(102) = {2, 0, 0, cl3};
Point(103) = {4, 0, 0, cl3};
Point(104) = {6, 0, 0, cl3};
Point(105) = {8, 0, 0, cl3};
Point(106) = {10, 0, 0, cl3};
Point(107) = {12, 0, 0, cl3};
Point(108) = {14, 0, 0, cl3};
Point(109) = {16, 0, 0, cl3};
Point(110) = {18, 0, 0, cl3};
Point(111) = {20, 0, 0, cl3};
Point(112) = {22, 0, 0, cl3};
Point(113) = {24, 0, 0, cl3};
Point(114) = {26, 0, 0, cl3};
Point(115) = {28, 0, 0, cl3};
Point(116) = {30, 0, 0, cl3};
Point(117) = {32, 0, 0, cl3};
Point(118) = {34, 0, 0, cl3};
Point(119) = {36, 0, 0, cl3};
Point(120) = {38, 0, 0, cl3};
Point(121) = {40, 0, 0, cl3};
Point(122) = {42, 0, 0, cl3};
Point(123) = {44, 0, 0, cl3};
Point(124) = {46, 0, 0, cl3};
Point(125) = {48, 0, 0, cl3};

// electrodes in north line
Point(201) = {0, 5, 0, cl3};
Point(202) = {2, 5, 0, cl3};
Point(203) = {4, 5, 0, cl3};
Point(204) = {6, 5, 0, cl3};
Point(205) = {8, 5, 0, cl3};
Point(206) = {10, 5, 0, cl3};
Point(207) = {12, 5, 0, cl3};
Point(208) = {14, 5, 0, cl3};
Point(209) = {16, 5, 0, cl3};
Point(210) = {18, 5, 0, cl3};
Point(211) = {20, 5, 0, cl3};
Point(212) = {22, 5, 0, cl3};
Point(213) = {24, 5, 0, cl3};
Point(214) = {26, 5, 0, cl3};
Point(215) = {28, 5, 0, cl3};
Point(216) = {30, 5, 0, cl3};
Point(217) = {32, 5, 0, cl3};
Point(218) = {34, 5, 0, cl3};
Point(219) = {36, 5, 0, cl3};
Point(220) = {38, 5, 0, cl3};
Point(221) = {40, 5, 0, cl3};
Point(222) = {42, 5, 0, cl3};
Point(223) = {44, 5, 0, cl3};
Point(224) = {46, 5, 0, cl3};
Point(225) = {48, 5, 0, cl3};

// electrodes in south line
Point(301) = {0, -5, 0, cl3};
Point(302) = {2, -5, 0, cl3};
Point(303) = {4, -5, 0, cl3};
Point(304) = {6, -5, 0, cl3};
Point(305) = {8, -5, 0, cl3};
Point(306) = {10, -5, 0, cl3};
Point(307) = {12, -5, 0, cl3};
Point(308) = {14, -5, 0, cl3};
Point(309) = {16, -5, 0, cl3};
Point(310) = {18, -5, 0, cl3};
Point(311) = {20, -5, 0, cl3};
Point(312) = {22, -5, 0, cl3};
Point(313) = {24, -5, 0, cl3};
Point(314) = {26, -5, 0, cl3};
Point(315) = {28, -5, 0, cl3};
Point(316) = {30, -5, 0, cl3};
Point(317) = {32, -5, 0, cl3};
Point(318) = {34, -5, 0, cl3};
Point(319) = {36, -5, 0, cl3};
Point(320) = {38, -5, 0, cl3};
Point(321) = {40, -5, 0, cl3};
Point(322) = {42, -5, 0, cl3};
Point(323) = {44, -5, 0, cl3};
Point(324) = {46, -5, 0, cl3};
Point(325) = {48, -5, 0, cl3};



//base of outer zone
Point(35) = {-125, 150, 0, cl1};
Point(36) = {175, 150, 0, cl1};
Point(37) = {175, -150, 0, cl1};
Point(38) = {-125, -150, 0, cl1};
Line(29) = {38, 35};
Line(30) = {35, 36};
Line(31) = {36, 37};
Line(32) = {37, 38};

//base of inner zone
Point(39) = {-25, 25, 0, cl2};
Point(40) = {75, 25, 0, cl2};
Point(41) = {75, -25, 0, cl2};
Point(42) = {-25, -25, 0, cl2};
Line(33) = {42, 39};
Line(34) = {39, 40};
Line(35) = {40, 41};
Line(36) = {41, 42};

//set base of inner zone
Extrude {0, 0, -20} {
  Line{33, 34, 35, 36};
}

//set base of outer zone
Extrude {0, 0, -100} {
  Line{29, 30, 31, 32};
}

Line Loop(69) = {53, 57, 61, 65};
Plane Surface(70) = {69};
Line Loop(71) = {29, 30, 31, 32};
Line Loop(72) = {33, 34, 35, 36};
Plane Surface(73) = {71, 72};
Line Loop(74) = {41, 45, 49, 37};
Plane Surface(75) = {74};
Plane Surface(76) = {72};
Surface Loop(77) = {70, 56, 73, 60, 64, 68, 40, 44, 48, 52, 75};

//define points in inner zone but out of circle
Point{101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 
201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 
301, 302, 303, 304, 305, 306, 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325
} In Surface{76};



Volume(78) = {77};
//Plane Surface(79) = {72};
//Surface Loop(80) = {40, 44, 48, 52, 75, 76};
//Surface Loop(81) = {6, 19, 23, 27};
//Volume(82) = {80, 81};


Surface Loop(79) = {76, 75, 44, 48, 52, 40};
Volume(80) = {79};
