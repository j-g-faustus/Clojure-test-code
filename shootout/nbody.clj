
(ns shootout.nbody
  (:require 
    [clojure.contrib.str-utils2 :as su] )
  (:gen-class) )

(comment"
Jovian planets N-body simulation.
From http://shootout.alioth.debian.org/u32q/benchmark.php?test=nbody&lang=all

Timing on MacBook Pro, N=50000000, java -Xms256M -Xmx256M
  real	0m50.777s
  user	0m51.494s

  $ java -version
  java version 1.6.0_20
  Java(TM) SE Runtime Environment (build 1.6.0_20-b02-279-10M3065)
  Java HotSpot(TM) 64-Bit Server VM (build 16.3-b01-279, mixed mode)

Timing for Java version - same machine, same N, same flags
    Java from http://shootout.alioth.debian.org/u64/benchmark.php?test=nbody&lang=java&sort=
  real	0m13.118s
  user	0m13.126s

See http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
and http://www.bestinclass.dk/index.clj/2010/03/functional-fluid-dynamics-in-clojure.html
for more on numeric performance tuning in Clojure.
")

(def *solar-mass* (* 4 Math/PI Math/PI))
(def *days-year* 365.24)

; Data for initial state
(def *data*
  ; names and tags for readability
  (list
    [:sun
     :x 0.0
     :y 0.0
     :z 0.0
     :vx 0.0
     :vy 0.0
     :vz 0.0
     :mass 1.0]

    [:jupiter
     :x 4.84143144246472090e+00
     :y -1.16032004402742839e+00
     :z -1.03622044471123109e-01
     :vx 1.66007664274403694e-03
     :vy 7.69901118419740425e-03
     :vz -6.90460016972063023e-05
     :mass 9.54791938424326609e-04]

    [:saturn
     :x 8.34336671824457987e+00
     :y 4.12479856412430479e+00
     :z -4.03523417114321381e-01
     :vx -2.76742510726862411e-03
     :vy 4.99852801234917238e-03
     :vz 2.30417297573763929e-05
     :mass 2.85885980666130812e-04]

    [:uranus
     :x 1.28943695621391310e+01
     :y -1.51111514016986312e+01
     :z -2.23307578892655734e-01
     :vx 2.96460137564761618e-03
     :vy 2.37847173959480950e-03
     :vz -2.96589568540237556e-05
     :mass 4.36624404335156298e-05]

    [:neptune
     :x 1.53796971148509165e+01
     :y -2.59193146099879641e+01
     :z 1.79258772950371181e-01
     :vx 2.68067772490389322e-03
     :vy 1.62824170038242295e-03
     :vz -9.51592254519715870e-05
     :mass 5.15138902046611451e-05] ))

; Convenience macros for array access
(defmacro x [body]
  `(aget ~body 0))
(defmacro y [body]
  `(aget ~body 1))
(defmacro z [body]
  `(aget ~body 2))
(defmacro vx [body]
  `(aget ~body 3))
(defmacro vy [body]
  `(aget ~body 4))
(defmacro vz [body]
  `(aget ~body 5))
(defmacro mass [body]
  `(aget ~body 6))

(defmacro x= [body val]
  `(aset ~body 0 ~val))
(defmacro y= [body val]
  `(aset ~body 1 ~val))
(defmacro z= [body val]
  `(aset ~body 2 ~val))
(defmacro vx= [body val]
  `(aset ~body 3 ~val))
(defmacro vy= [body val]
  `(aset ~body 4 ~val))
(defmacro vz= [body val]
  `(aset ~body 5 ~val))
(defmacro mass= [body val]
  `(aset ~body 6 ~val))

(defmacro x+= [body val]
  (let [b# body]
    `(x= ~b# (+ (x ~b#) ~val))))
(defmacro y+= [body val]
  (let [b# body]
    `(y= ~b# (+ (y ~b#) ~val))))
(defmacro z+= [body val]
  (let [b# body]
    `(z= ~b# (+ (z ~b#) ~val))))
(defmacro vx+= [body val]
  (let [b# body]
    `(vx= ~b# (+ (vx ~b#) ~val))))
(defmacro vy+= [body val]
  (let [b# body]
    `(vy= ~b# (+ (vy ~b#) ~val))))
(defmacro vz+= [body val]
  (let [b# body]
    `(vz= ~b# (+ (vz ~b#) ~val))))

; Convert to Java Object[] of double[]
(defn convert-data []
  ; to Object[]
  (into-array Object
    ; to double[]
    (map (fn [body]
           (let [^doubles b (into-array Double/TYPE
                     (map last (partition 2 (rest body)))) ]
             (vx= b (* (vx b) *days-year*))
             (vy= b (* (vy b) *days-year*))
             (vz= b (* (vz b) *days-year*))
             (mass= b (* (mass b) *solar-mass*))
             b))
      *data*)))


(def *bodies* (convert-data))

(defn reset-state []
  (def *bodies* (convert-data)))

(defn init-state []
  "Initializes state"
  (reset-state)
  (let [[px py pz] (reduce (fn [[px py pz] ^doubles b]
                             [(+ px (* (vx b) (mass b)))
                              (+ py (* (vy b) (mass b)))
                              (+ pz (* (vz b) (mass b))) ])
                     [0.0 0.0 0.0]
                     *bodies*)
        sun (aget *bodies* 0)
        mass (mass sun) ]
    (vx= sun (/ (- px) mass))
    (vy= sun (/ (- py) mass))
    (vz= sun (/ (- pz) mass))))

(defn energy
  "Total energy for current state"
  ([] ; start
    (energy *bodies*))
  ([bodies] ; sum
    (if-not bodies 0.0
      (+ (energy (first bodies) (next bodies))
        (energy (next bodies)))))
  ([body nbodies] ; one body
    (let [v-sq (+ (* (vx body) (vx body))
                 (* (vy body) (vy body))
                 (* (vz body) (vz body)))
          e (* 1/2 (mass body) v-sq)  ; kinetic energy: 1/2 * m * v^2
          f (fn [e nbody]
              (let [dx (- (x body) (x nbody))
                    dy (- (y body) (y nbody))
                    dz (- (z body) (z nbody))
                    dist (Math/sqrt
                           (+ (* dx dx) (* dy dy) (* dz dz))) ]
                (- e (/ (* (mass body) (mass nbody))
                       dist)))) ]
      (reduce f e nbodies))))

(defn advance 
  "Move system one dt timestep forwards"
  [dt]
  (let [#^objects bodies *bodies* ; x2 improvement
        len (int (alength bodies))
        dt (double dt) ]  ; x10 improvement 
    ; update velocity
    (loop [i (int 0)]
      (when (< i len)
        (let [body (doubles (aget bodies i))
              bx (x body)
              by (y body)
              bz (z body) ]
          (loop [j (int (unchecked-inc i)) ]
            (when (< j len)
              (let [nbody (doubles (aget bodies j))
                    dx (- bx (x nbody))
                    dy (- by (y nbody))
                    dz (- bz (z nbody))
                    dsq (+ (* dx dx)
                          (+ (* dy dy)
                            (* dz dz)))
                    mag (/ dt
                          (* dsq (Math/sqrt dsq))) ]
                (let [mult (- (* (mass nbody) mag)) ] ; subtraction
                  (vx+= body (* dx mult))
                  (vy+= body (* dy mult))
                  (vz+= body (* dz mult)))
                (let [mult (* (mass body) mag) ]
                  (vx+= nbody (* dx mult))
                  (vy+= nbody (* dy mult))
                  (vz+= nbody (* dz mult)))
                (recur (unchecked-inc j)) )))

          (recur (unchecked-inc i)) )))

    ; update position
    (loop [i (int 0)]
      (when (< i len)
        (let [body (doubles (aget bodies i)) ]
          (x+= body (* dt (vx body)))
          (y+= body (* dt (vy body)))
          (z+= body (* dt (vz body))))
        (recur (unchecked-inc i))))))

(defn -main [& args]
  (let [[n x] args]
    (when (= x "-wait")
      (println "Press enter to start")
      (read-line))
    (init-state)
    (when (not= x :bench)
      (println (format "%.9f" (energy))))
    (dotimes [_ (int (Integer/parseInt n))]
      (advance 0.01))
    (when (not= x :bench)
      (println (format "%.9f" (energy))))))

; Current timing n=5,000,000: 5.5 - 5.8 s, avg. around 5.6 after 3-4 runs.
; Update: 5.3-5.4 s, after typecasting inner loop j to int
; Update: 5.2-5.4 s, after skipping double typecast in aget!
; Update: 5.1-5.3 s, after ??? - seems to be related to appropriate level of
; typecasts. (Not casting all doubles, only some of them.)
(defn time-test []
  (time (-main "1000000")))

(defn time-full []
  (time (-main "50000000")))


;(defn bench-test []
; (cc/bench (-main "1000000" :bench)))

; Verify correct output for n = 1000
; Data: http://shootout.alioth.debian.org/u32q/iofile.php?test=nbody&file=output
(defn- verify []
  (let [res (with-out-str (-main "1000"))
        [e-init e-end] ["-0.169075164" "-0.169087605"]
        [init end] (su/split res #"\s+") ]
    (print "Start state ")
    (if (= init e-init)
      (println "OK")
      (println "failed:" init "should be" e-init))
    (print "End state ")
    (if (= end e-end)
      (println "OK")
      (println "failed:" end "should be:" e-end))))
