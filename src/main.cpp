#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.rfind("}]");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
  double result = 0.0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * pow(x, i);
  }
  return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
  assert(xvals.size() == yvals.size());
  assert(order >= 1 && order <= xvals.size() - 1);
  Eigen::MatrixXd A(xvals.size(), order + 1);

  for (int i = 0; i < xvals.size(); i++) {
    A(i, 0) = 1.0;
  }

  for (int j = 0; j < xvals.size(); j++) {
    for (int i = 0; i < order; i++) {
      A(j, i + 1) = A(j, i) * xvals(j);
    }
  }

  auto Q = A.householderQr();
  auto result = Q.solve(yvals);
  return result;
}

//transform to car x coordinate function
vector<double> transform_x(vector<double> ptsx,
				vector<double> ptsy,
				double px,
				double py,
				double psi) {
   double dx;
   double dy;
   vector<double> tx;
   
   //Transform ptsx to car x coordinate
   for (int j = 0; j < ptsx.size(); j++) {
	dx = ptsx[j] - px;
	dy = ptsy[j] - py;
	
	tx.push_back(dx * cos(psi) + dy * sin(psi));
	
   }
   return tx; 
}

//transform to car y coordinate function
vector<double> transform_y(vector<double> ptsx,
				vector<double> ptsy,
				double px,
				double py,
				double psi) {
   double dx;
   double dy;
   vector<double> ty;
   
   //Transform ptsy to car y coordinate
   for (int k = 0; k < ptsy.size(); k++) {
	dx = ptsx[k] - px;
	dy = ptsy[k] - py;
	ty.push_back(-dx * sin(psi) + dy * cos(psi));
	
   }
   return ty;
}
int main() {
  uWS::Hub h;

  // MPC is initialized here!
  MPC mpc;

  h.onMessage([&mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    string sdata = string(data).substr(0, length);
    cout << sdata << endl;
    if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
      string s = hasData(sdata);
      if (s != "") {
        auto j = json::parse(s);
        string event = j[0].get<string>();
        if (event == "telemetry") {
          // j[1] is the data JSON object
          vector<double> ptsx = j[1]["ptsx"];
          vector<double> ptsy = j[1]["ptsy"];
          double px = j[1]["x"];
          double py = j[1]["y"];
          double psi = j[1]["psi"];
          double v = j[1]["speed"];

          /*
          * TODO: Calculate steering angle and throttle using MPC.
          *
          * Both are in between [-1, 1].
          *
          */
          double steer_value;
          double throttle_value;

	  //Transform to car coordinates
	  vector<double> ptsx_co = transform_x(ptsx, ptsy, px, py, psi);
	  vector<double> ptsy_co = transform_y(ptsx, ptsy, px, py, psi);

	  //convert pts data array into Eigen vector array
	  Eigen::VectorXd ptsxv(ptsx_co.size());
	  Eigen::VectorXd ptsyv(ptsy_co.size());
	  ptsxv = Eigen::VectorXd::Map(ptsx_co.data(), ptsx_co.size());
	  ptsyv = Eigen::VectorXd::Map(ptsy_co.data(), ptsy_co.size());

	  auto coeffs = polyfit(ptsxv, ptsyv, 3);

	  //since x = 0, y = 0, and angle = 0
          //cte is calculated by evaluating at x(-1) and subtracting y.
	  double cte = polyeval(coeffs, 0);

	  //check derivation is correct since angle = psi = 0
	  double epsi = -atan(coeffs[1]);
          std::cout << "cte:  " << cte << "	v:  " << v << "	epsi:  " << epsi << std::endl;

	  Eigen::VectorXd state(6);
	  double next_x;
	  double next_y;
	  double next_psi;
	  double next_v;
	  double next_cte;
	  double next_epsi;
	  const double Lf = 2.67;

	  //predicted state is 100ms ahead.
	  //The time delay for 100ms is 0.1s.
	  double delta_t = 0.1;
	  double prev_delta = mpc.prev_delta;
	  double prev_a = mpc.prev_a;

	  //if x=0, y = 0, and psi = 0, the global kinetmatic model next state variables are:
	  next_x = v * delta_t;
	  next_y = 0;
	  next_psi = -v * (prev_delta/Lf) * delta_t;
	  next_v = v + prev_a * delta_t;
	  next_cte = cte + v * CppAD::sin(epsi) * delta_t;
	  next_epsi = epsi + v * (prev_delta/Lf) * delta_t;

  	  state << next_x, next_y, next_psi, next_v, next_cte, next_epsi;
	  
          //The MPC uses an optimizer(Ipopt solver) to find the control inputs that minimize
          //the cost function.
    	  auto vars = mpc.Solve(state, coeffs);

    	  steer_value = vars[0]/(deg2rad(25) * Lf);
	  throttle_value = vars[1];
	  std::cout << "steer_value: " << steer_value << "	throttle_value: " << throttle_value << std::endl;

	  //keep the previous values for later use
	  mpc.prev_delta = steer_value;
	  mpc.prev_a = throttle_value;

          json msgJson;
          // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
          // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
          msgJson["steering_angle"] = steer_value;
          msgJson["throttle"] = throttle_value;

          //Display the MPC predicted trajectory(line projection) in green
          vector<double> mpc_x_vals;
          vector<double> mpc_y_vals;

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Green line
	 for (int j = 2; j < vars.size(); j++) {
            if(j % 2 == 0){
            	mpc_x_vals.push_back(vars[j]);
            } else {
              	mpc_y_vals.push_back(vars[j]);
            }
	    //std::cout << "vars[" << j << "]: " << vars[j] << std::endl;
          }
          
          msgJson["mpc_x"] = mpc_x_vals;
          msgJson["mpc_y"] = mpc_y_vals;

          //Display the waypoints/reference line (polynomial fitted reference line) in yellow
          vector<double> next_x_vals;
          vector<double> next_y_vals;
          
          next_x_vals.resize(ptsxv.size());
          next_y_vals.resize(ptsyv.size());
          

          //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
          // the points in the simulator are connected by a Yellow line

          for (int k = 0; k < ptsxv.size(); k++) {
            next_x_vals[k] = ptsxv[k];
            next_y_vals[k] = ptsyv[k];
          }
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"steer\"," + msgJson.dump() + "]";
          std::cout << msg << std::endl;
          // Latency
          // The purpose is to mimic real driving conditions where
          // the car does actuate the commands instantly.
          //
          // Feel free to play around with this value but should be to drive
          // around the track with 100ms latency.
          //
          // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
          // SUBMITTING.
          this_thread::sleep_for(chrono::milliseconds(100));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
