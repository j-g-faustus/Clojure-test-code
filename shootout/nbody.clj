
(ns shootout.nbody
  (:require 
    [clojure.contrib.str-utils2 :as su] )
  (:gen-class)
  )

(comment
  Jovian planets N-body simulation.
  From http://shootout.alioth.debian.org/u32q/benchmark.php?test=nbody&lang=all
  )

(def *solar-mass* (* 4 Math/PI Math/PI))
(def *days-year* 365.24)

; Data for initial state
(def *data*
  ; list constants with names and tags for readability
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

; Convert to Java Object[] of double[]
(defn convert-data []
  ; to Object[]
  (into-array Object
    ; to double[]
    (map (fn [[x y z vx vy vz mass]]
           (into-array Double/TYPE
             (list x y z
               (* vx *days-year*)
               (* vy *days-year*)
               (* vz *days-year*)
               (* mass *solar-mass*))))
      ; to just numbers
      (map (fn [planet]
             (map last (partition 2 (rest planet))))
        *data*))))


(def *bodies* (convert-data))

(defn reset []
  (def *bodies* (convert-data)))

; Macros for speed.
; See http://clj-me.cgrand.net/2009/10/15/multidim-arrays/
; and http://www.bestinclass.dk/index.clj/2010/03/functional-fluid-dynamics-in-clojure.html
(defmacro aget!
  ([array i]
    `(double (aget ~(vary-meta array assoc :tag 'doubles) (int ~i)))))

(defmacro aset!
  ([array i v]
    `(aset ~(vary-meta array assoc :tag 'doubles) (int ~i) (double ~v))))

; Convenience macros for array access
(defmacro x [body]
  `(aget! ~body 0))
(defmacro y [body]
  `(aget! ~body 1))
(defmacro z [body]
  `(aget! ~body 2))
(defmacro vx [body]
  `(aget! ~body 3))
(defmacro vy [body]
  `(aget! ~body 4))
(defmacro vz [body]
  `(aget! ~body 5))
(defmacro mass [body]
  `(aget! ~body 6))

(defmacro x= [body val]
  `(aset! ~body (int 0) ~val))
(defmacro y= [body val]
  `(aset! ~body (int 1) ~val))
(defmacro z= [body val]
  `(aset! ~body (int 2) ~val))
(defmacro vx= [body val]
  `(aset! ~body (int 3) ~val))
(defmacro vy= [body val]
  `(aset! ~body (int 4) ~val))
(defmacro vz= [body val]
  `(aset! ~body (int 5) ~val))

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


(defn init 
  "Initializes state"
  []
  (let [#^objects bodies *bodies*
        len (int (alength bodies)) ]
    (loop [i (int 0)
           px (double 0.0)
           py (double 0.0)
           pz (double 0.0) ]
      (if (>= i len)
        (let [sun (doubles (aget bodies 0)) ]
          (vx= sun (/ (- px) (double *solar-mass*)))
          (vy= sun (/ (- py) (double *solar-mass*)))
          (vz= sun (/ (- pz) (double *solar-mass*))))
        (let [body (doubles (aget bodies i))
              bmass (double (mass body)) ]
          (recur
            (int (inc i))
            (+ px (* (vx body) bmass))
            (+ py (* (vy body) bmass))
            (+ pz (* (vz body) bmass))))))))

(defn advance 
  "Move system one dt timestep forwards"
  [dt]
  (let [#^objects bodies *bodies*
        len (int (alength bodies))
        dt (double dt) ]
      ; update velocities
      (loop [i (int 0)]
        (when (< i len)
          (let [body (doubles (aget bodies i))
                bx (double (x body))
                by (double (y body))
                bz (double (z body)) ]
            (loop [j (unchecked-inc i)]
              (when (< j len)
                (let [nbody (doubles (aget bodies j))
                      dx (double (- bx (x nbody)))
                      dy (double (- by (y nbody)))
                      dz (double (- bz (z nbody)))
                      dsq (double
                            (+ (* dx dx)
                              (+ (* dy dy)
                                (* dz dz))))
                      dist (double (Math/sqrt dsq))
                      mag (double (/ dt (* dsq dist))) ]
                  (let [mult (double (* (mass nbody) mag))]
                    (vx+= body (- (* dx mult)))
                    (vy+= body (- (* dy mult)))
                    (vz+= body (- (* dz mult))))
                  (let [mult (double (* (mass body) mag))]
                    (vx+= nbody (* dx mult))
                    (vy+= nbody (* dy mult))
                    (vz+= nbody (* dz mult))))
                (recur (unchecked-inc j)))))

          (recur (unchecked-inc i))))

       ; update position
      (loop [i (int 0)]
        (when (< i len)
          (let [body (doubles (aget bodies i)) ]
            (x+= body (* dt (vx body)))
            (y+= body (* dt (vy body)))
            (z+= body (* dt (vz body))))
          (recur (unchecked-inc i))))))


(defn energy 
  "Returns total energy for current state"
  []
  (let [#^objects bodies *bodies*
        len (int (alength bodies)) ]
    (loop [i (int 0)
           e (double 0.0)]
      (if-not (< i len) e
      (let [body (doubles (aget bodies i))
            ne (double
                 (+ e
                   (* (* 0.5 (mass body))
                     (+ (* (vx body)  (vx body) )
                       (+ (* (vy body) (vy body))
                         (* (vz body) (vz body)))))))
            nne (double
                  (loop [j (int (inc i))
                         nne ne ]
                    (if-not (< j len) nne
                      (let [nbody (doubles (aget bodies j))
                            dx (double (- (x body) (x nbody)))
                            dy (double (- (y body) (y nbody)))
                            dz (double (- (z body) (z nbody)))
                            dist (double
                                   (Math/sqrt
                                     (+ (* dx dx)
                                       (+ (* dy dy)
                                         (* dz dz))))) ]
                        (recur (int (inc j))
                          (- nne
                            (/ (* (mass body) (mass nbody))
                              dist))))))) ]
        (recur (int (inc i)) nne))))))



(defn -main
  "Run simulation"
  [& args]
  (let [n (Integer/parseInt (first args)) ]
    (reset)
    (init)
    (println (format "%.9f" (energy)))
    (dotimes [i (int n)]
      (advance 0.01))
    (println (format "%.9f" (energy)))))

; Verify correct output for n = 1000
; Data: http://shootout.alioth.debian.org/u32q/iofile.php?test=nbody&file=output
(defn- verify []
  (let [res (with-out-str (-main "1000"))
       [init end] (su/split res #"\s+") ]
    (if (= init "-0.169075164")
      (println "start state: OK")
      (println "start state failed:" init))
    (if (= end "-0.169087605")
      (println "end state: OK")
      (println "end state failed:" end "should be:" "-0.169087605"))))
