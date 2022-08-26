#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#define PI 3.14159265
#define WIN -100000000
#define MUL 10000
#define LL long long

using namespace std;

/**
 * Auto-generated code below aims at helping you parse
 * the standard input according to the problem statement.
 **/

class Mt {
  int angles[40];
  int anglesc[40];

public:
  void init() {
    for (int a = -90; a <= 90; a += 5) {
      angles[(90 + a) / 5] = sin(PI * a / 180) * MUL;
      anglesc[(90 + a) / 5] = cos(PI * a / 180) * MUL;
    }
  }
  int get(int a) {
    // cerr<<"sin("<<a<<")="<<angles[(90+a)/5]<<endl;
    return angles[(90 + a) / 5];
  }
  int getc(int a) { return anglesc[(90 + a) / 5]; }
};

struct pt {
  int x, y;
};

inline int area(pt a, pt b, pt c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

inline bool intersect_1(int a, int b, int c, int d) {
  if (a > b)
    swap(a, b);
  if (c > d)
    swap(c, d);
  return max(a, c) <= min(b, d);
}

bool intersect(pt a, pt b, pt c, pt d) {
  return intersect_1(a.x, b.x, c.x, d.x) && intersect_1(a.y, b.y, c.y, d.y) &&
         area(a, b, c) * area(a, b, d) <= 0 &&
         area(c, d, a) * area(c, d, b) <= 0;
}

class Shuttle {
private:
  LL X;
  LL Y;
  int hSpeed;
  int vSpeed;
  int fuel;
  int rotate;
  int power;

public:
  void init(LL p1, LL p2, int p3, int p4, int p5, int p6, int p7) {
    X = p1;
    Y = p2;
    hSpeed = p3;
    vSpeed = p4;
    fuel = p5;
    rotate = p6;
    power = p7;
    // cerr<<"init Y:"<<Y<<endl;
  }

  int getAngle() { return rotate; }
  LL getX() { return X; }
  LL getY() { return Y; }
  int getVX() { return hSpeed; }
  int getVY() { return vSpeed; }
  bool Possible() { return abs(hSpeed) <= 2000 && abs(vSpeed) <= 4000; }

  Shuttle next(int n_ang, int ipw, double gM, Mt m_sin) {
    Shuttle rez;
    double pw =
        (abs(ipw - power) > 1 ? (ipw < power ? power - 1 : power + 1) : ipw);
    int d_ang = min(15, abs(rotate - n_ang));
    int ang = (rotate > n_ang ? rotate - d_ang : rotate + d_ang);
    LL nx = (X + hSpeed + (int)(((pw * 100 * (m_sin.get(-ang)) / MUL) / 2)));
    int vx = (hSpeed + ((pw * 100 * m_sin.get(-ang)) / MUL));
    // cerr<<"vx:"<<hSpeed<<"+(("<<pw<<"*100*"<<m_sin.get(ang)<<")/"<<MUL<<")) =
    // "<<vx<<endl;
    LL ny =
        (Y + vSpeed + (int)(((pw * 100 * m_sin.getc(-ang)) / MUL) - gM) / 2);
    int vy = vSpeed + (int)(((pw * 100 * m_sin.getc(-ang)) / MUL) - gM);

    rez.init(nx, ny, vx, vy, fuel - ipw, ang, ipw);
    return rez;
  }
};

class Delta {
public:
  int dx, dy, dvx, dvy;
  void init(int _dx, int _dy, int _dvx, int _dvy) {
    dx = _dx;
    dy = _dy;
    dvx = _dvx;
    dvy = _dvy;
  }
};

class StDelta {
public:
  Delta mas[4 * 40 + (90 + 90) / 5];
  void add(int a, int b, int dx, int dy, int dvx, int dvy) {
    mas[b * 37 + (90 + a) / 5].init(dx, dy, dvx, dvy);
  }
  Delta get(int a, int b) { return mas[b * 37 + (90 + a) / 5]; }
};

class Leaf {
public:
  int parent_id;
  Delta n;
  int angle;
  int power;
};

class Graph {
private:
  const int MaxVx = 1800;
  const int MaxVy = 3800;
  const double gMars = 371.1;
  int stepping;
  StDelta sa_sd;
  vector<int> sx, sy;
  int surfaceN;
  int ml1, ml2, my;
  int old_y = -10000000, old_x = 0;
  int it;
  Shuttle nomad;
  Shuttle test;
  int rez_a, rez_p, c_a, c_p, max_cost;
  int m_c15;
  Mt m_sin;
  LL alt[7000];
  int best_way;
  LL cc1, cc2, cc3;

public:
  void init(int sN) {
    surfaceN = sN;
    sx.resize(sN);
    sy.resize(sN);
    it = 0;
    m_sin.init();
  }

  int hLand(){
    return my;
  }
  void add_point(int landX, int landY) {
    landX *= 100;
    landY *= 100;
    if (it > 0) {
      for (int xx = sx[it - 1] / 100; xx < landX / 100; xx++) {
        alt[xx] = sy[it - 1] +
                  (xx * 100 - sx[it - 1]) * (landY - sy[it - 1]) /
                      (landX - sx[it - 1]) +
                  500;
      }
    }
    sx[it] = landX;
    sy[it] = landY;
    it++;
    if (old_x >= 0) {
      if (landY == old_y) {
        if (landX - old_x > 99999) {
          ml1 = old_x;
          ml2 = landX;
          my = landY;
          old_x = -1;
          cerr << "land is " << ml1 << ":" << ml2 << "," << my << endl;
        } else {
          old_y = landY;
          old_x = landX;
        }
      } else {
        old_y = landY;
        old_x = landX;
      }
    }
  }

  LL get_cost(Shuttle n, bool fl_dbg=false) {
    LL _mid_x = (ml1 + ml2) / 2;

    LL _x = n.getX();
    LL _y = n.getY();
    LL _vx = n.getVX();
    LL _vy = n.getVY();
    LL avx = abs(_vx);
    LL avy = abs(_vy);
    LL cost0 = 0;
    if (_x <= ml2 - 8500 && _x >= ml1 + 8500 && _y <= my && avx <= 1900 &&
        avy <= 3500 && n.getAngle() == 0)
      cost0 = WIN;
    LL cost1 =
        sqrt((_x + 5*_vx - _mid_x) * (_x  + 5*_vx - _mid_x) + (_y  + 5*_vy - my) * (_y  + 5*_vy - my)) *10;
    if ((_x < ml1 && _vx <= 0) || (_x > ml2 && _vx >= 0))
      avx += 10000;
    LL cost2 = (avx > 2000 ? (avx > 2500 ? 300 * avx : 40 * avx) : avx);
    LL cost3 = (avy > 3500 ? (avy > 4000 ? 1500 * avy : 200 * avy) : avy);
    cc1 = cost1;
    cc2 = cost2;
    cc3 = cost3;
    LL cost4 = abs(n.getAngle())*2;
    if(_y<my)cost4+=10000;
    return cost0 + cost1 + cost2 + cost3 + cost4;
  }

  LL get_cost(Leaf n, bool fl_dbg=false) {
    LL _mid_x = (ml1 + ml2) / 2;

    LL _x = n.n.dx;
    LL _y = n.n.dy;
    LL _vx = n.n.dvx;
    LL _vy = n.n.dvy;
    LL avx = abs(_vx);
    LL avy = abs(_vy);
    LL cost0 = 0;

    if (_x <= ml2 - 8500 && _x >= ml1 + 8500 && _y <= my && avx <= 1900 &&
        avy <= 3500 && n.angle == 0)
      cost0 = WIN;
    LL cost1 =
        sqrt((_x + 5*_vx - _mid_x) * (_x  + 5*_vx - _mid_x) + (_y  + 5*_vy - my) * (_y  + 5*_vy - my)) *10;
    if ((_x < ml1 && _vx <= 0) || (_x > ml2 && _vx >= 0))
      avx += 10000;
    LL cost2 = (avx > 2000 ? (avx > 2500 ? 300 * avx : 40 * avx) : avx);
    LL cost3 = (avy > 3500 ? (avy > 4000 ? 1500 * avy : 200 * avy) : avy);
    cc1 = cost1;
    cc2 = cost2;
    cc3 = cost3;
    LL cost4 = abs(n.angle)*2;
    if(_y<my)cost4+=10000;
    return cost0 + cost1 + cost2 + cost3 + cost4;
  }

  LL found_way(vector<Leaf> *g, int st_id, int step) {
    int sz = (*g).size();
    LL best_cost = -WIN*5;
    best_way = st_id;
    for (int j = st_id; j < sz; j++) {
      Leaf l = (*g)[j];
      int power = l.power;
      int rotate = l.angle;
      Delta n = l.n;
      int old_p = -2342;
      for (int p = max(0, power - 1); p <= min(4, power + 1); p++)
        for (int a = max(-90, rotate - 15);
             a <= min(90, rotate + 15); a += 5) {
          stepping++;
          Delta da = sa_sd.get(a, p);
          Delta na;
          Leaf nl;
          nl.parent_id = j;
          nl.angle = a;
          nl.power = p;
          na.dx = n.dx;
          na.dy = n.dy;
          na.dvx = n.dvx;
          na.dvy = n.dvy;
          for (int time = 0; time < (a==rotate ? 20 : 3); time++) {
            na.dx += na.dvx + da.dx;
            na.dy += na.dvy + da.dy;
            na.dvx += da.dvx;
            na.dvy += da.dvy;
          }
          nl.n = na;
          if (na.dx <= ml2 - 8500 && na.dx >= ml1 + 8500 &&
              na.dy <= my + 15000 && abs(na.dvx) <= 1900 &&
              abs(na.dvy) <= 3900 && a == 0) {
            best_cost = WIN;
            best_way = j;
            return WIN;
          } else if (na.dx > 0 && na.dx < 700000 && na.dy < 300000) {
            if (na.dy > alt[na.dx / 100] || (na.dx <= ml2 - 8500 && na.dx >= ml1 + 8500))
                if(step==3){
                    LL t_cost = get_cost(nl);
                    if(t_cost<best_cost){
                        best_cost = t_cost;
                        best_way = j;
                    }
                }
                else {
                    (*g).push_back(nl);
                }
            else {
                if(-WIN*2<best_cost){
                    best_cost = -WIN*2;
                    best_way = j;
                }
            }
            } else {
                if(-WIN<best_cost){
                    best_cost = -WIN;
                    best_way = j;
                }
            }
        }
    }
    if(step==3){
      int c_id = best_way;
      while (c_id!=-1){
        c_id = (*g)[c_id].parent_id;
      }
      return best_cost;
    }
    best_cost = found_way(g, sz, step + 1);
    if(step>1)best_way = (*g)[best_way].parent_id;
    return best_cost;
  }

  void first(LL X, LL Y, int hSpeed, int vSpeed, int fuel, int rotate,
             int power) {
    vector<Leaf> graph;
    Shuttle t;
    for (int p = 0; p <= 4; p++) {
      for (int a = -90; a <= 90; a += 5) {
        t.init(0, 0, 0, 0, 0, a, p);
        t = t.next(a, p, gMars, m_sin);
        sa_sd.add(a, p, t.getX(), t.getY(), t.getVX(), t.getVY());
      }
    }
    nomad.init(X, Y, hSpeed, vSpeed, fuel, rotate, power);
    rez_a = rotate;
    rez_p = power;
    int id = -1;
    Leaf l;
    Delta inp;
    inp.init(X, Y, hSpeed, vSpeed);
    l.n = inp;
    l.parent_id = id;
    l.angle = rotate;
    l.power = power;
    graph.push_back(l);
    stepping = 0;
    LL t_cost = found_way(&graph, 0, 1);
    rez_a = graph[best_way].angle;
    rez_p = graph[best_way].power;
    if(t_cost==WIN){
      if(abs(vSpeed)>3000)rez_p = 4;
      if(abs(Y+vSpeed-my)<12000)rez_a=0;
    }
    //compare(X, Y, hSpeed, vSpeed, fuel, rotate, power);
  }

  void out() {
    nomad = nomad.next(rez_a, rez_p, gMars, m_sin);
    cout << rez_a << " " << rez_p << endl;
  }
};

int main() {
  int surfaceN; // the number of points used to draw the surface of Mars.
  cin >> surfaceN;
  cin.ignore();
  Graph sl;
  sl.init(surfaceN);

  for (int i = 0; i < surfaceN; i++) {
    int landX; // X coordinate of a surface point. (0 to 6999)
    int landY; // Y coordinate of a surface point. By linking all the points
               // together in a sequential fashion, you form the surface of
               // Mars.
    cin >> landX >> landY;
    cin.ignore();
    sl.add_point(landX, landY);
  }
  bool fl_init = false;
  int turn=0;
  // game loop
  while (1) {
    int X;
    int Y;
    int hSpeed; // the horizontal speed (in m/s), can be negative.
    int vSpeed; // the vertical speed (in m/s), can be negative.
    int fuel;   // the quantity of remaining fuel in liters.
    int rotate; // the rotation angle in degrees (-90 to 90).
    int power;  // the thrust power (0 to 4).
    cin >> X >> Y >> hSpeed >> vSpeed >> fuel >> rotate >> power;
    cin.ignore();
    if (!fl_init) {
      sl.first(X * 100, Y * 100, hSpeed * 100, vSpeed * 100, fuel, rotate,
               power);
      // sl.compare(X*100,Y*100,hSpeed*100,vSpeed*100,fuel,rotate,power);
      fl_init = true;
    } else {
      if(sl.hLand()>16000&&turn++<80){
        cout<<"10 4"<<endl;
        continue;
      }
      else sl.first(X * 100, Y * 100, hSpeed * 100, vSpeed * 100, fuel, rotate,
                 power);
    }

    // Write an action using cout. DON'T FORGET THE "<< endl"
    // To debug: cerr << "Debug messages..." << endl;

    // rotate power. rotate is the desired rotation angle. power is the desired
    // thrust power.
    sl.out();
  }
}
