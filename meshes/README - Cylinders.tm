<TeXmacs|1.0.7.19>

<style|<tuple|generic|maxima>>

<\body>
  <section|Tuning the extrusion for the cylinders>

  We'd like to refine the mesh as we approach the midpoint of the cylinder's
  long axis. For this we we fix the distance between layers of the extrusion
  to follow a simple law. Assume the cylinder's long axis is given by the
  segment <math|<around*|[|-a,a|]>> and that we want <math|N> subdivisions at
  each side of the middle point <math|x=0>, at points with <math|x>
  coordinate <math|x<rsub|k>>. We may use the following recursive formula\ 

  <\equation>
    x<rsub|k+1>=x<rsub|k>+d<around*|(|x<rsub|k>|)><label|eq:def-seq>
  </equation>

  where

  <\equation>
    d<around*|(|x<rsub|k>|)>=<frac|a-<big|sum><rsub|j=1><rsup|k-1>d<around*|(|x<rsub|j>|)>|N-k+1>*<frac|N|a>=<frac|-x<rsub|k>|N-k+1>*<frac|N|a><label|eq:def-dist-seq>
  </equation>

  is the distance between <math|x<rsub|k+1>> and <math|x<rsub|k>>. If we set
  <math|x<rsub|1>=-a>, then this choice of <math|d> ensures that the next
  term is

  <\equation*>
    x<rsub|2>=x<rsub|1>+1.
  </equation*>

  Using <eqref|eq:def-seq> and <eqref|eq:def-dist-seq> we have the simple
  expression

  <\equation*>
    x<rsub|k+1>=x<rsub|k>*<around*|(|1-<frac|N|a*<around*|(|N-k+1|)>>|)>=c<around*|(|k|)>*x<rsub|k>.
  </equation*>

  Because of the way the <tt|Layers> command works in <name|Gmsh>, we need to
  calculate the normalized width of each layer:

  <\equation*>
    <wide|x|\<wide-bar\>><rsub|k>=<frac|a-x<rsub|k>|2a>.
  </equation*>

  Thus, the <tt|Extrude> command will look like:

  <\indent>
    <\code>
      Extrude {

      \ \ Surface (surface_id_base);

      \ \ Layers { {1, <with|color|dark green|... n times total ...>, 1},

      \ \ \ \ \ \ \ \ \ \ \ { <with|color|dark
      green|<math|<wide|x|\<wide-bar\>><rsub|1>>>, <with|color|dark
      green|<math|<wide|x|\<wide-bar\>><rsub|2>>>, ..., <with|color|dark
      green|<math|<wide|x|\<wide-bar\>><rsub|n>>> }};

      \ \ Recombine;

      };
    </code>
  </indent>

  <\session|scheme|default>
    <\input|Scheme] >
      (define N 20) ; \ The number of subdivisions at each side of 0
    </input>

    <\input|Scheme] >
      (define a 5) ; The half-width of the cylinder
    </input>

    <\input|Scheme] >
      (define (factor k)

      \ \ "A helper function to calculate the factor c(k)"

      \ \ (if (== k (+ N 1)) 0 ; safeguard

      \ \ \ \ \ \ (- 1 (/ N (* a (+ N 1 (- 0 k)))))))
    </input>

    <\input|Scheme] >
      (define (term k)

      \ \ "Returns the k'th term x_k"

      \ \ (if (== k 1)\ 

      \ \ \ \ \ \ (- 0 a) \ ; First term must be -a to have first distance ==
      1

      \ \ \ \ \ \ (* (term (- k 1)) (factor (- k 1)))))
    </input>

    <\input|Scheme] >
      (define (range from to)

      \ \ "A useful function"

      \ \ (if (== from to) (list to)

      \ \ (cons from (range ((if (\<less\> from to) + -) from 1) to))))))
    </input>

    <\input|Scheme] >
      (define seq1 (map term (range 1 20)))
    </input>

    <\unfolded-io|Scheme] >
      seq1
    <|unfolded-io>
      (-5 -4 -60/19 -140/57 -1820/969 -455/323 -1001/969 -715/969 -165/323
      -110/323 -70/323 -42/323 -70/969 -35/969 -5/323 -5/969 -1/969 0 0 0)
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  Hmm... I don't understand those zeroes at the end, maybe it's a limitation
  of the integer arithmetic in this <scheme> implementation. For the moment
  I'll just truncate the result:

  <\session|scheme|default>
    <\input|Scheme] >
      (define seq1 (cDDDr seq1)) ; HACK
    </input>

    <\unfolded-io|Scheme] >
      seq1
    <|unfolded-io>
      (-5 -4 -60/19 -140/57 -1820/969 -455/323 -1001/969 -715/969 -165/323
      -110/323 -70/323 -42/323 -70/969 -35/969 -5/323 -5/969 -1/969)
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  We now mirror this sequence for the positive interval
  <math|<around*|[|0,a|]>>:

  <\session|scheme|default>
    <\input|Scheme] >
      (define seq2 (map (lambda (x) (* (- 0 1) x)) (reverse seq1)))
    </input>

    <\input|Scheme] >
      (define seq (append seq1 seq2))
    </input>

    <\unfolded-io|Scheme] >
      seq
    <|unfolded-io>
      (-5 -4 -60/19 -140/57 -1820/969 -455/323 -1001/969 -715/969 -165/323
      -110/323 -70/323 -42/323 -70/969 -35/969 -5/323 -5/969 -1/969 1/969
      5/969 5/323 35/969 70/969 42/323 70/323 110/323 165/323 715/969
      1001/969 455/323 1820/969 140/57 60/19 4 5)
    </unfolded-io>

    <\unfolded-io|Scheme] >
      (length seq)
    <|unfolded-io>
      34
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  The list <scm|seq> contains the points <math|x<rsub|k>,k=1,\<ldots\>,2*N-1>
  where the layers must be created, so we still need to calculate the
  normalized <math|x<rsub|k>>. First a test:

  <\session|scheme|default>
    <\input|Scheme] >
      (define (dx l r)

      \ \ "Compute the differences between elements of a list of numbers"

      \ \ (if (== 1 (length l)) r

      \ \ \ \ \ \ (cons (- (cadr l) (car l)) (dx (cdr l) r))))
    </input>

    <\unfolded-io|Scheme] >
      (dx seq '())
    <|unfolded-io>
      (1 16/19 40/57 560/969 455/969 364/969 286/969 220/969 55/323 40/323
      28/323 56/969 35/969 20/969 10/969 4/969 2/969 4/969 10/969 20/969
      35/969 56/969 28/323 40/323 55/323 220/969 286/969 364/969 455/969
      560/969 40/57 16/19 1)
    </unfolded-io>

    <\unfolded-io|Scheme] >
      (apply + (dx seq '()))
    <|unfolded-io>
      10
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  The sum is correctly equal to the interval's width <math|2*a>. Now we
  compute the normalized sequence

  <\equation*>
    <wide|x|\<wide-bar\>><rsub|k>=<frac|a-x<rsub|k>|2a>.
  </equation*>

  <\session|scheme|default>
    <\input|Scheme] >
      (define (normalize x) (/ (+ a x) (* 2 a))) seq))
    </input>

    <\unfolded-io|Scheme] >
      (map normalize seq)
    <|unfolded-io>
      (0 1/10 7/38 29/114 605/1938 116/323 1922/4845 413/969 145/323 301/646
      309/646 1573/3230 955/1938 481/969 161/323 484/969 2422/4845 2423/4845
      485/969 162/323 488/969 983/1938 1657/3230 337/646 345/646 178/323
      556/969 2923/4845 207/323 1333/1938 85/114 31/38 9/10 1)
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  We may try to plot this inside <TeXmacs>:

  <center|<with|gr-mode|<tuple|edit|point>|gr-frame|<tuple|scale|0.840795cm|<tuple|0.5gw|0.160031gh>>|gr-geometry|<tuple|geometry|0.613349par|0.160018par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|1|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-grid-aspect-props|<tuple|<tuple|axes|#808080>|<tuple|1|#c0c0c0>|<tuple|10|#e0e0ff>>|gr-grid-aspect|<tuple|<tuple|axes|#808080>|<tuple|1|#c0c0c0>>|magnify|0.840896413466164|<graphics||<point|-5.0|0.0>|<point|-4.0|0.1>|<point|-3.15789473684211|0.184210526315789>|<point|-2.45614035087719|0.254385964912281>|<point|-1.87822497420021|0.312177502579979>|<point|-1.40866873065015|0.359133126934984>|<point|-1.03302373581011|0.396697626418989>|<point|-0.737874097007224|0.426212590299278>|<point|-0.510835913312694|0.448916408668731>|<point|-0.340557275541796|0.46594427244582>|<point|-0.21671826625387|0.478328173374613>|<point|-0.130030959752322|0.486996904024768>|<point|-0.0722394220846233|0.492776057791538>|<point|-0.0361197110423117|0.496388028895769>|<point|-0.0154798761609907|0.498452012383901>|<point|-0.00515995872033024|0.499484004127967>|<point|-0.00103199174406605|0.499896800825593>|<point|0.00103199174406605|0.500103199174407>|<point|0.00515995872033024|0.500515995872033>|<point|0.0154798761609907|0.501547987616099>|<point|0.0361197110423117|0.503611971104231>|<point|0.0722394220846233|0.507223942208462>|<point|0.130030959752322|0.513003095975232>|<point|0.21671826625387|0.521671826625387>|<point|0.340557275541796|0.53405572755418>|<point|0.510835913312694|0.551083591331269>|<point|0.737874097007224|0.573787409700722>|<point|1.03302373581011|0.603302373581011>|<point|1.40866873065015|0.640866873065015>|<point|1.87822497420021|0.687822497420021>|<point|2.45614035087719|0.745614035087719>|<point|3.15789473684211|0.81578947368421>|<point|4.0|0.9>|<point|5.0|1.0>|<line|<point|-5.0|0.0>|<point|-4.0|0.1>|<point|-3.15789473684211|0.184210526315789>|<point|-2.45614035087719|0.254385964912281>|<point|-1.87822497420021|0.312177502579979>|<point|-1.40866873065015|0.359133126934984>|<point|-1.03302373581011|0.396697626418989>|<point|-0.737874097007224|0.426212590299278>|<point|-0.510835913312694|0.448916408668731>|<point|-0.340557275541796|0.46594427244582>|<point|-0.21671826625387|0.478328173374613>|<point|-0.130030959752322|0.486996904024768>|<point|-0.0722394220846233|0.492776057791538>|<point|-0.0361197110423117|0.496388028895769>|<point|-0.0154798761609907|0.498452012383901>|<point|-0.00515995872033024|0.499484004127967>|<point|-0.00103199174406605|0.499896800825593>|<point|0.00103199174406605|0.500103199174407>|<point|0.00515995872033024|0.500515995872033>|<point|0.0154798761609907|0.501547987616099>|<point|0.0361197110423117|0.503611971104231>|<point|0.0722394220846233|0.507223942208462>|<point|0.130030959752322|0.513003095975232>|<point|0.21671826625387|0.521671826625387>|<point|0.340557275541796|0.53405572755418>|<point|0.510835913312694|0.551083591331269>|<point|0.737874097007224|0.573787409700722>|<point|1.03302373581011|0.603302373581011>|<point|1.40866873065015|0.640866873065015>|<point|1.87822497420021|0.687822497420021>|<point|2.45614035087719|0.745614035087719>|<point|3.15789473684211|0.81578947368421>|<point|4.0|0.9>|<point|5.0|1.0>>>>>

  Interestingly, the sequence obtained is a cubic function of <math|k>. Here
  we display the position of the subdivisions

  <center|<with|gr-mode|<tuple|edit|point>|gr-frame|<tuple|scale|0.706795cm|<tuple|252699tmpt|0.508435gh>>|gr-geometry|<tuple|geometry|0.546685par|0.533336par|center>|gr-grid|<tuple|cartesian|<point|0|0>|1>|gr-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-aspect|<tuple|<tuple|axes|none>|<tuple|1|none>|<tuple|1|none>>|gr-edit-grid|<tuple|cartesian|<point|0|0>|1>|gr-edit-grid-old|<tuple|cartesian|<point|0|0>|1>|gr-grid-aspect-props|<tuple|<tuple|axes|#808080>|<tuple|1|#c0c0c0>|<tuple|10|#e0e0ff>>|gr-grid-aspect|<tuple|<tuple|axes|#808080>|<tuple|1|#c0c0c0>>|magnify|0.707106777750335|<graphics||<point|-5.0|-5.0>|<point|-4.70588235294118|-4.0>|<point|-4.41176470588235|-3.15789473684211>|<point|-4.11764705882353|-2.45614035087719>|<point|-3.82352941176471|-1.87822497420021>|<point|-3.52941176470588|-1.40866873065015>|<point|-3.23529411764706|-1.03302373581011>|<point|-2.94117647058824|-0.737874097007224>|<point|-2.64705882352941|-0.510835913312694>|<point|-2.35294117647059|-0.340557275541796>|<point|-2.05882352941176|-0.21671826625387>|<point|-1.76470588235294|-0.130030959752322>|<point|-1.47058823529412|-0.0722394220846233>|<point|-1.17647058823529|-0.0361197110423117>|<point|-0.882352941176471|-0.0154798761609907>|<point|-0.588235294117647|-0.00515995872033024>|<point|-0.294117647058824|-0.00103199174406605>|<point|0.0|0.00103199174406605>|<point|0.294117647058824|0.00515995872033024>|<point|0.588235294117647|0.0154798761609907>|<point|0.882352941176471|0.0361197110423117>|<point|1.17647058823529|0.0722394220846233>|<point|1.47058823529412|0.130030959752322>|<point|1.76470588235294|0.21671826625387>|<point|2.05882352941176|0.340557275541796>|<point|2.35294117647059|0.510835913312694>|<point|2.64705882352941|0.737874097007224>|<point|2.94117647058824|1.03302373581011>|<point|3.23529411764706|1.40866873065015>|<point|3.52941176470588|1.87822497420021>|<point|3.82352941176471|2.45614035087719>|<point|4.11764705882353|3.15789473684211>|<point|4.41176470588235|4.0>|<point|4.70588235294118|5.0>|<line|<point|-5.0|-5.0>|<point|-4.70588235294118|-4.0>|<point|-4.41176470588235|-3.15789473684211>|<point|-4.11764705882353|-2.45614035087719>|<point|-3.82352941176471|-1.87822497420021>|<point|-3.52941176470588|-1.40866873065015>|<point|-3.23529411764706|-1.03302373581011>|<point|-2.94117647058824|-0.737874097007224>|<point|-2.64705882352941|-0.510835913312694>|<point|-2.35294117647059|-0.340557275541796>|<point|-2.05882352941176|-0.21671826625387>|<point|-1.76470588235294|-0.130030959752322>|<point|-1.47058823529412|-0.0722394220846233>|<point|-1.17647058823529|-0.0361197110423117>|<point|-0.882352941176471|-0.0154798761609907>|<point|-0.588235294117647|-0.00515995872033024>|<point|-0.294117647058824|-0.00103199174406605>|<point|0.0|0.00103199174406605>|<point|0.294117647058824|0.00515995872033024>|<point|0.588235294117647|0.0154798761609907>|<point|0.882352941176471|0.0361197110423117>|<point|1.17647058823529|0.0722394220846233>|<point|1.47058823529412|0.130030959752322>|<point|1.76470588235294|0.21671826625387>|<point|2.05882352941176|0.340557275541796>|<point|2.35294117647059|0.510835913312694>|<point|2.64705882352941|0.737874097007224>|<point|2.94117647058824|1.03302373581011>|<point|3.23529411764706|1.40866873065015>|<point|3.52941176470588|1.87822497420021>|<point|3.82352941176471|2.45614035087719>|<point|4.11764705882353|3.15789473684211>|<point|4.41176470588235|4.0>|<point|4.70588235294118|5.0>>>>>

  <\session|scheme|default>
    <\input|Scheme] >
      (define (point x y)

      \ \ `(point ,(number-\<gtr\>string (exact-\<gtr\>inexact x))\ 

      \ \ \ \ \ \ \ \ \ \ ,(number-\<gtr\>string (exact-\<gtr\>inexact y))))
    </input>

    <\input|Scheme] >
      (define points (map point seq (map normalize seq)))
    </input>

    <\input|Scheme] >
      (define points2

      \ \ (with c (/ (length seq) 2)

      \ \ \ \ (with l (map (lambda (x) (* (/ x 17) 5)) (range (- 0 c) (- c
      1)))

      \ \ \ \ \ \ (map point l seq))))
    </input>

    <\input|Scheme] >
      (define (update-plots)

      \ \ "MEGA HACK"

      \ \ (let* ((grs (select (buffer-tree) '(:* graphics)))

      \ \ \ \ \ \ \ \ \ (gr1 (car grs))

      \ \ \ \ \ \ \ \ \ (gr2 (cadr grs)))

      \ \ \ \ (tree-set gr1 (stree-\<gtr\>tree `(graphics ""\ 

      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ,@points\ 

      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (line
      ,@points))))

      \ \ \ \ (tree-set gr2 (stree-\<gtr\>tree `(graphics ""\ 

      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ ,@points2\ 

      \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ \ (line
      ,@points2)))))

      \ \ #t)
    </input>

    <\unfolded-io|Scheme] >
      (update-plots)
    <|unfolded-io>
      #t
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  Finally, we transform this into a string ready to paste into the
  <name|Gmsh> file:

  <\session|scheme|default>
    <\unfolded-io|Scheme] >
      (string-concatenate\ 

      \ (list-intersperse (map number-\<gtr\>string (map normalize seq)) ",
      "))
    <|unfolded-io>
      "0, 1/10, 7/38, 29/114, 605/1938, 116/323, 1922/4845, 413/969, 145/323,
      301/646, 309/646, 1573/3230, 955/1938, 481/969, 161/323, 484/969,
      2422/4845, 2423/4845, 485/969, 162/323, 488/969, 983/1938, 1657/3230,
      337/646, 345/646, 178/323, 556/969, 2923/4845, 207/323, 1333/1938,
      85/114, 31/38, 9/10, 1"
    </unfolded-io>

    <\input|Scheme] >
      \;
    </input>
  </session>

  The following <name|Gmsh> code will create the cylinder.

  <\indent>
    <\code>
      Point(1) = {0,0,-5,0.5};

      \;

      Point(2) = {1,0,-5,0.5};

      Point(3) = {0,1,-5,0.5};

      Point(4) = {-1,0,-5,0.5};

      Point(5) = {0,-1,-5,0.5};

      \;

      Circle(1) = {2,1,3};

      Circle(2) = {3,1,4};

      Circle(3) = {4,1,5};

      Circle(4) = {5,1,2};

      \;

      Line Loop(5) = {1,2,3,4};

      Plane Surface(6) = {5};

      Recombine Surface{6};

      Physical Surface("bottom") = {6};

      \;

      cyl[] = Extrude {0,0,10} {

      \ \ Surface{6};

      \ \ Layers { { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
      },

      \ \ \ \ \ \ \ \ \ \ \ {1/10, 7/38, 29/114, 605/1938, 116/323,
      1922/4845, 413/969, 145/323, 301/646, 309/646, 1573/3230, 955/1938,
      481/969, 161/323, 484/969, 2422/4845, 2423/4845, 485/969, 162/323,
      488/969, 983/1938, 1657/3230, 337/646, 345/646, 178/323, 556/969,
      2923/4845, 207/323, 1333/1938, 85/114, 31/38, 9/10, 1 }};

      \ \ Recombine;

      };
    </code>
  </indent>

  And this is the result after some visualization options have been set:

  <center|<image|cylinder master.eps|10cm|||>>
</body>

<\initial>
  <\collection>
    <associate|info-flag|detailed>
    <associate|preamble|false>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|?>>
    <associate|eq:def-dist-seq|<tuple|2|?>>
    <associate|eq:def-seq|<tuple|1|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Tuning
      the extrusion for the cylinders> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>