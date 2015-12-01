
// plot two curves together with legends
function plot2(points0, points1) {

    var plot_options = {
        legend:{
            backgroundOpacity: 0.1,
            noColumns: 0,
            position: "ne"
        },
        yaxes: [{ min:-0.05, max:1.19}],
        xaxes: [{ min:-0.8, max:12}],
    };

    $.plot('#div-plot', [
        {data: points0, color: "blue", label: 'CF 0'},
        {data: points1, color: "red", label: 'CF 1'}
    ], plot_options);

}

(function(angular) {
    'use strict';
    angular.module('ode', []).controller('odeController', function(){
        // Initial values
        this.freeze0 = false;
        this.isPandemic = false;

        this.PreCF = 1;
        this.Event=0.5;
        this.ER=0.5;
        this.SC=0.5;
        this.PR=0.5;
        this.PM=0.5;
        this.PVID=0.5;
        this.CFdplt=3;
        this.Edecay=4;
        this.ERflow=1;
        this.PRflow=1;
        this.SCflow=1;

        // Settings
        this.varlist = ["CF", "Event", "ER", "SC", "PR", "PreCF", "PVID", "PM"];


        this.getW = function getW0(){
            var W = [];
            for(var i=0; i < this.varlist.length; i++){
                var item = this.varlist[i];
                var cmd = "W[" + i + "] = parseFloat(this." + item +");";
                eval(cmd);
            }
            return W;
        };

        this.update = function update(){
            this.CF = this.PreCF;
            var W0 = this.getW();
            var sol = this.solve(W0, 0, 12);
            this.rslt1 = this.calculate_resistance(sol);
            // Plot
            var points1 = [];
            points1.push([-0.9, this.PreCF]);
            for (var i = 0; i < sol.CF.length; i++){
                points1.push([sol.tspan[i], sol.CF[i]]);
            }
            if (! this.freeze0){
                this.points0 = points1;
                this.rslt0 = this.rslt1;
            }
            plot2(this.points0, points1);


        };


        this.fdW = function fdW(t, W){
            var local = [];
            for(var i=0; i < this.varlist.length; i++){
                var item = this.varlist[i];
                local[item] = parseFloat(W[i]);
            }

            // Replenish
            var CF_drop = Math.max( local.PreCF - local.CF, 0 );
            var SC_flow_rate = this.SCflow * local.CF * local.SC * CF_drop;
            var PR_flow_rate = this.PRflow * local.PR * CF_drop;
            var ER_flow_rate = this.ERflow * local.ER * CF_drop;
            var CF_replenish_rate = SC_flow_rate + PR_flow_rate + ER_flow_rate;


            // CF depletion rate
            var Event_decay_rate = this.Edecay * local.Event * (local.PM + local.PVID)/2;
            if (this.isPandemic){
                var mu = 2.0; // event peak (mean)
                var sigma = 1.0; // event spread (std)
                var pi = Math.PI;
                Event_decay_rate *= (1/Math.sqrt(2*pi*sigma*sigma)) * Math.exp(-(t-mu)*(t-mu) / (2*sigma*sigma));
            }

            var CF_depletion_rate = local.CF * Event_decay_rate;

            // Derivatives
            var d = [];
            d.SC = -SC_flow_rate;
            d.PR = - PR_flow_rate;
            d.ER = - ER_flow_rate;
            d.Event = - Math.max(Event_decay_rate, 0);
            d.CF = CF_replenish_rate - this.CFdplt * CF_depletion_rate;
            d.PreCF = 0;
            d.PVID = 0;
            d.PM = 0;

            var dW = [];
            for(var i=0; i < this.varlist.length; i++){
                var cmd = "dW[" + i + "] = d." + this.varlist[i] + ";";
                eval(cmd);
            }
            return dW;
        };


        this.solve = function solve(W0, tmin, tmax){
            // 4th order Runge-Kutta method
            var sol = {};
            sol.CF = [];
            sol.tspan = [];

            var W = W0;
            var t = tmin;
            while(t <= tmax){
                sol.CF.push(W[0]);
                sol.tspan.push(t);
                var dt = (t<1)? 0.05 : 0.2;
                W = this.int_one_step(W, t, dt);
                t += dt;
            }
            return sol;
        };

        this.int_one_step = function int_one_step(W, t, dt){
            function aXbY(a, X, b, Y){
                // Calculate a * X + b * Y
                var Z = [];
                for(var i=0; i < X.length; i++){
                    Z[i] = a * X[i] + b * Y[i];
                }
                return Z;
            }

            var k1 = this.fdW(t, W);
            var k2 = this.fdW(t+0.5*dt, aXbY(1, W, 0.5*dt, k1));
            var k3 = this.fdW(t+0.5*dt, aXbY(1, W, 0.5*dt, k2));
            var k4 = this.fdW(t+dt, aXbY(1, W, dt, k3));
            var k12 = aXbY(1, k1, 2, k2);
            var k34 = aXbY(2, k3, 1, k4);
            var k1234 = aXbY(1, k12, 1, k34);
            var W = aXbY(1, W, dt/6, k1234);
            return W;
        };

        this.calculate_resistance = function calculate_resistance(sol){
            var rslt = {};
            var CFmin = d3.min(sol.CF);
            var Idxmin = sol.CF.indexOf(CFmin);
            var CF0 = sol.CF[0];
            rslt.Resistance = 100 * CFmin / CF0;
            // Find T half using bisection
            var CFhalf = 0.5 * (CFmin + CF0);
            var Idx1 = Idxmin;
            var Idx2 = sol.CF.length - 1;
            while(Idx2 - Idx1 > 1){
                var Idx0 = Math.round( 0.5 * ( Idx1 + Idx2) );
                if(sol.CF[Idx0] == CFhalf){
                    break;
                }else{
                    if(sol.CF[Idx0] > CFhalf){
                        Idx2 = Idx0;
                    }else{
                        Idx1 = Idx0;
                    }
                }
            }
            rslt.tHalf = sol.tspan[Idx0];
            rslt.Resilience = d3.mean(sol.CF);
            rslt.Recovery = 1 / rslt.tHalf;

            return rslt;
        };


        this.update();

    });
})(window.angular);