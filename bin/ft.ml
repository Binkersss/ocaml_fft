(*  ISC License
 *
 *  Copyright (c) 2026 Nathaniel Chappelle
 *
 *  Permission to use, copy, modify, and/or distribute this software for any
 *  purpose with or without fee is hereby granted, provided that the above
 *  copyright notice and this permission notice appear in all copies.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH
 *  REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY
 *  AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT,
 *  INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 *  LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR
 *  OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR
 *  PERFORMANCE OF THIS SOFTWARE.
 *
 *  NAS Parallel Benchmarks FT (3D FFT PDE solver), sequential OCaml port.
 *  Build:  ocamlopt -unsafe -inline 100 -I +unix unix.cmxa ft.ml \
 *            -cclib -lnuma -o ft.native
 *  Run:    FT_FAR=0 ./ft.native [S|W|A|B|C]
 *)

(* ---------- problem classes ---------- *)
type params = { nx : int; ny : int; nz : int; niter : int; name : string }

let params_of = function
  | "S" -> { nx = 64; ny = 64; nz = 64; niter = 6; name = "S" }
  | "W" -> { nx = 128; ny = 128; nz = 32; niter = 6; name = "W" }
  | "A" -> { nx = 256; ny = 256; nz = 128; niter = 6; name = "A" }
  | "B" -> { nx = 512; ny = 256; nz = 256; niter = 20; name = "B" }
  | "C" -> { nx = 512; ny = 512; nz = 512; niter = 20; name = "C" }
  | c -> failwith ("unknown class " ^ c)

(* ---------- constants ---------- *)
let fftblock = 16
let pad = 18 (* FFTBLOCKPAD: row stride of the scratch blocks *)
let seed = 314159265.0
let a_mult = 1220703125.0
let alpha = 1.0e-6
let pi = 3.141592653589793238

(* ---------- NPB linear congruential RNG (x_{k+1} = a*x_k mod 2^46) ---------- *)
let r23 = 0.5 ** 23.0
let r46 = r23 *. r23
let t23 = 2.0 ** 23.0
let t46 = t23 *. t23
let trunc x = float_of_int (int_of_float x)

let randlc (x : float ref) a =
  let a1 = trunc (r23 *. a) in
  let a2 = a -. (t23 *. a1) in
  let x1 = trunc (r23 *. !x) in
  let x2 = !x -. (t23 *. x1) in
  let t1 = (a1 *. x2) +. (a2 *. x1) in
  let t2 = trunc (r23 *. t1) in
  let z = t1 -. (t23 *. t2) in
  let t3 = (t23 *. z) +. (a2 *. x2) in
  let t4 = trunc (r46 *. t3) in
  x := t3 -. (t46 *. t4);
  r46 *. !x

(* fills y.(0..n-1), updates seed *)
let vranlc n (x_seed : float ref) a (y : float array) =
  let a1 = trunc (r23 *. a) in
  let a2 = a -. (t23 *. a1) in
  let x = ref !x_seed in
  for i = 0 to n - 1 do
    let x1 = trunc (r23 *. !x) in
    let x2 = !x -. (t23 *. x1) in
    let t1 = (a1 *. x2) +. (a2 *. x1) in
    let t2 = trunc (r23 *. t1) in
    let z = t1 -. (t23 *. t2) in
    let t3 = (t23 *. z) +. (a2 *. x2) in
    let t4 = trunc (r46 *. t3) in
    x := t3 -. (t46 *. t4);
    y.(i) <- r46 *. !x
  done;
  x_seed := !x

(* a^exponent mod 2^46 by repeated squaring through randlc *)
let ipow46 a exponent =
  if exponent = 0 then 1.0
  else begin
    let q = ref a and r = ref 1.0 and n = ref exponent in
    while !n > 1 do
      let n2 = !n / 2 in
      if n2 * 2 = !n then (
        ignore (randlc q !q);
        n := n2)
      else (
        ignore (randlc r !q);
        n := !n - 1)
    done;
    ignore (randlc r !q);
    !r
  end

(* ---------- complex 3D array: separate re/im float arrays ----------
   index of element (k,j,i) = (k*ny + j)*nx + i   (k: z, j: y, i: x)      *)
type c3 = { re : float array; im : float array }

let make_c3 ~far n =
  { re = Array.make ~far n 0.0; im = Array.make ~far n 0.0 }

let gen_initial_conditions ~far g =
  let plane = g.nx * g.ny in
  let u = make_c3 ~far (plane * g.nz) in
  let tmp = Array.make (2 * plane) 0.0 in
  let start = ref seed in
  (* the original's first jump is a^0 = 1, which leaves [start] unchanged *)
  let an = ipow46 a_mult (2 * plane) in
  for k = 0 to g.nz - 1 do
    let x0 = ref !start in
    vranlc (2 * plane) x0 a_mult tmp;
    for p = 0 to plane - 1 do
      u.re.((k * plane) + p) <- tmp.(2 * p);
      u.im.((k * plane) + p) <- tmp.((2 * p) + 1)
    done;
    ignore (randlc start an)
  done;
  u

(* ---------- helpers ---------- *)
let ilog2 n =
  let rec go n acc = if n <= 1 then acc else go (n / 2) (acc + 1) in
  go n 0

(* roots of unity table: U[0] = (m, 0); level j occupies U[2^j .. 2^(j+1)-1] *)
let fft_init n ur ui =
  let m = ilog2 n in
  ur.(0) <- float_of_int m;
  ui.(0) <- 0.0;
  let ku = ref 1 and ln = ref 1 in
  for _ = 0 to m - 1 do
    let t = pi /. float_of_int !ln in
    for i = 0 to !ln - 1 do
      let ti = float_of_int i *. t in
      ur.(i + !ku) <- cos ti;
      ui.(i + !ku) <- sin ti
    done;
    ku := !ku + !ln;
    ln := 2 * !ln
  done

(* ---------- Stockham FFT on a (n x pad) block, columns 0..ny-1 ---------- *)
let fftz2 is l m n ny ur ui xr xi yr yi =
  let n1 = n / 2 in
  let lk = 1 lsl (l - 1) in
  let li = 1 lsl (m - l) in
  let lj = 2 * lk in
  for i = 0 to li - 1 do
    let i11 = i * lk in
    let i12 = i11 + n1 in
    let i21 = i * lj in
    let i22 = i21 + lk in
    let u1r = ur.(li + i) in
    let u1i = if is >= 1 then ui.(li + i) else -.ui.(li + i) in
    for k = 0 to lk - 1 do
      for j = 0 to ny - 1 do
        let a = ((i11 + k) * pad) + j and b = ((i12 + k) * pad) + j in
        let c = ((i21 + k) * pad) + j and d = ((i22 + k) * pad) + j in
        let x11r = xr.(a) and x11i = xi.(a) in
        let x21r = xr.(b) and x21i = xi.(b) in
        yr.(c) <- x11r +. x21r;
        yi.(c) <- x11i +. x21i;
        let dr = x11r -. x21r and di = x11i -. x21i in
        yr.(d) <- (u1r *. dr) -. (u1i *. di);
        yi.(d) <- (u1r *. di) +. (u1i *. dr)
      done
    done
  done

(* result ends up in x *)
let cfftz ur ui is m n xr xi yr yi =
  let l = ref 1 in
  while !l <= m do
    fftz2 is !l m n fftblock ur ui xr xi yr yi;
    if !l + 1 <= m then fftz2 is (!l + 1) m n fftblock ur ui yr yi xr xi;
    l := !l + 2
  done;
  if m mod 2 = 1 then
    for j = 0 to n - 1 do
      for i = 0 to fftblock - 1 do
        xr.((j * pad) + i) <- yr.((j * pad) + i);
        xi.((j * pad) + i) <- yi.((j * pad) + i)
      done
    done

(* FFT along x (dimension of size nx), 16 y-rows at a time *)
let cffts1 g ur ui is (xin : c3) (xout : c3) =
  let nx, ny, nz = (g.nx, g.ny, g.nz) in
  let logn = ilog2 nx in
  let y0r = Array.make (nx * pad) 0.0 and y0i = Array.make (nx * pad) 0.0 in
  let y1r = Array.make (nx * pad) 0.0 and y1i = Array.make (nx * pad) 0.0 in
  for k = 0 to nz - 1 do
    let jj = ref 0 in
    while !jj <= ny - fftblock do
      for j = 0 to fftblock - 1 do
        for i = 0 to nx - 1 do
          let s = (((k * ny) + j + !jj) * nx) + i in
          y0r.((i * pad) + j) <- xin.re.(s);
          y0i.((i * pad) + j) <- xin.im.(s)
        done
      done;
      cfftz ur ui is logn nx y0r y0i y1r y1i;
      for j = 0 to fftblock - 1 do
        for i = 0 to nx - 1 do
          let s = (((k * ny) + j + !jj) * nx) + i in
          xout.re.(s) <- y0r.((i * pad) + j);
          xout.im.(s) <- y0i.((i * pad) + j)
        done
      done;
      jj := !jj + fftblock
    done
  done

(* FFT along y (size ny), 16 x-columns at a time *)
let cffts2 g ur ui is (xin : c3) (xout : c3) =
  let nx, ny, nz = (g.nx, g.ny, g.nz) in
  let logn = ilog2 ny in
  let m = max nx ny in
  let y0r = Array.make (m * pad) 0.0 and y0i = Array.make (m * pad) 0.0 in
  let y1r = Array.make (m * pad) 0.0 and y1i = Array.make (m * pad) 0.0 in
  for k = 0 to nz - 1 do
    let ii = ref 0 in
    while !ii <= nx - fftblock do
      for j = 0 to ny - 1 do
        for i = 0 to fftblock - 1 do
          let s = (((k * ny) + j) * nx) + i + !ii in
          y0r.((j * pad) + i) <- xin.re.(s);
          y0i.((j * pad) + i) <- xin.im.(s)
        done
      done;
      cfftz ur ui is logn ny y0r y0i y1r y1i;
      for j = 0 to ny - 1 do
        for i = 0 to fftblock - 1 do
          let s = (((k * ny) + j) * nx) + i + !ii in
          xout.re.(s) <- y0r.((j * pad) + i);
          xout.im.(s) <- y0i.((j * pad) + i)
        done
      done;
      ii := !ii + fftblock
    done
  done

(* FFT along z (size nz), 16 x-columns at a time *)
let cffts3 g ur ui is (xin : c3) (xout : c3) =
  let nx, ny, nz = (g.nx, g.ny, g.nz) in
  let logn = ilog2 nz in
  let m = max nx nz in
  let y0r = Array.make (m * pad) 0.0 and y0i = Array.make (m * pad) 0.0 in
  let y1r = Array.make (m * pad) 0.0 and y1i = Array.make (m * pad) 0.0 in
  for j = 0 to ny - 1 do
    let ii = ref 0 in
    while !ii <= nx - fftblock do
      for k = 0 to nz - 1 do
        for i = 0 to fftblock - 1 do
          let s = (((k * ny) + j) * nx) + i + !ii in
          y0r.((k * pad) + i) <- xin.re.(s);
          y0i.((k * pad) + i) <- xin.im.(s)
        done
      done;
      cfftz ur ui is logn nz y0r y0i y1r y1i;
      for k = 0 to nz - 1 do
        for i = 0 to fftblock - 1 do
          let s = (((k * ny) + j) * nx) + i + !ii in
          xout.re.(s) <- y0r.((k * pad) + i);
          xout.im.(s) <- y0i.((k * pad) + i)
        done
      done;
      ii := !ii + fftblock
    done
  done

(* forward: x, y, z -> output; inverse: z, y, x -> output *)
let fft g ur ui dir (input : c3) (output : c3) =
  if dir = 1 then begin
    cffts1 g ur ui 1 input input;
    cffts2 g ur ui 1 input input;
    cffts3 g ur ui 1 input output
  end
  else begin
    cffts3 g ur ui (-1) input input;
    cffts2 g ur ui (-1) input input;
    cffts1 g ur ui (-1) input output
  end

(* ---------- time evolution ---------- *)
let compute_indexmap ~far g =
  let { nx; ny; nz; _ } = g in
  let idx = Array.make ~far (nx * ny * nz) 0 in
  let fold i n = ((i + (n / 2)) mod n) - (n / 2) in
  for k = 0 to nz - 1 do
    let kk = fold k nz in
    for j = 0 to ny - 1 do
      let jj = fold j ny in
      for i = 0 to nx - 1 do
        let ii = fold i nx in
        idx.((((k * ny) + j) * nx) + i) <- (ii * ii) + (jj * jj) + (kk * kk)
      done
    done
  done;
  idx

let compute_exp ~far g =
  let expmax =
    g.niter * ((g.nx * g.nx / 4) + (g.ny * g.ny / 4) + (g.nz * g.nz / 4))
  in
  let ex = Array.make ~far (expmax + 1) 0.0 in
  ex.(0) <- 1.0;
  ex.(1) <- exp (-4.0 *. alpha *. pi *. pi);
  for i = 2 to expmax do
    ex.(i) <- ex.(i - 1) *. ex.(1)
  done;
  ex

let evolve (u0 : c3) (u1 : c3) t indexmap ex =
  for p = 0 to Array.length indexmap - 1 do
    let e = ex.(t * indexmap.(p)) in
    u1.re.(p) <- u0.re.(p) *. e;
    u1.im.(p) <- u0.im.(p) *. e
  done

let checksum g sums iter (u : c3) =
  let { nx; ny; nz; _ } = g in
  let cr = ref 0.0 and ci = ref 0.0 in
  for j = 1 to 1024 do
    let q = j mod nx and r = 3 * j mod ny and s = 5 * j mod nz in
    let p = (((s * ny) + r) * nx) + q in
    cr := !cr +. u.re.(p);
    ci := !ci +. u.im.(p)
  done;
  let ntotal = float_of_int (nx * ny * nz) in
  sums.(iter) <- (!cr /. ntotal, !ci /. ntotal)

(* ---------- verification (reference values from NPB, scaled by 1e2) ---------- *)
let ref_values = function
  | "S" ->
      ( [
          5.546087004964;
          5.546385409189;
          5.546148406171;
          5.545423607415;
          5.544255039624;
          5.542683411902;
        ],
        [
          4.845363331978;
          4.865304269511;
          4.883910722336;
          4.901273169046;
          4.917475857993;
          4.932597244941;
        ] )
  | "W" ->
      ( [
          5.673612178944;
          5.631436885271;
          5.594024089970;
          5.560698047020;
          5.530898991250;
          5.504159734538;
        ],
        [
          5.293246849175;
          5.282149986629;
          5.270996558037;
          5.260027904925;
          5.249400845633;
          5.239212247086;
        ] )
  | "A" ->
      ( [
          5.046735008193;
          5.059412319734;
          5.069376896287;
          5.077892868474;
          5.085233095391;
          5.091487099959;
        ],
        [
          5.114047905510;
          5.098809666433;
          5.098144042213;
          5.101336130759;
          5.104914655194;
          5.107917842803;
        ] )
  | "B" ->
      ( [
          5.177643571579;
          5.154521291263;
          5.146409228649;
          5.142378756213;
          5.139626667737;
          5.137423460082;
          5.135547056878;
          5.133910925466;
          5.132470705390;
          5.131197729984;
          5.130070319283;
          5.129070537032;
          5.128182883502;
          5.127393733383;
          5.126691062020;
          5.126064276004;
          5.125504076570;
          5.125002331720;
          5.124551951846;
          5.124146770029;
        ],
        [
          5.077803458597;
          5.088249431599;
          5.096208912659;
          5.101023387619;
          5.103976610617;
          5.105948019802;
          5.107404165783;
          5.108576573661;
          5.109577278523;
          5.110460304483;
          5.111252433800;
          5.111968077718;
          5.112616233064;
          5.113203605551;
          5.113735928093;
          5.114218460548;
          5.114656139760;
          5.115053595966;
          5.115415130407;
          5.115744692211;
        ] )
  | "C" ->
      ( [
          5.195078707457;
          5.155422171134;
          5.144678022222;
          5.140150594328;
          5.137550426810;
          5.135811056728;
          5.134569343165;
          5.133651975661;
          5.132955192805;
          5.132410471738;
          5.131971141679;
          5.131605205716;
          5.131290734194;
          5.131012720314;
          5.130760908195;
          5.130528295923;
          5.130310107773;
          5.130103090133;
          5.129905029333;
          5.129714421109;
        ],
        [
          5.149019699238;
          5.127578201997;
          5.122251847514;
          5.121090289018;
          5.121143685824;
          5.121496764568;
          5.121870921893;
          5.122193250322;
          5.122454735794;
          5.122663649603;
          5.122830879827;
          5.122965869718;
          5.123075927445;
          5.123166486553;
          5.123241541685;
          5.123304037599;
          5.123356167976;
          5.123399592211;
          5.123435588985;
          5.123465164008;
        ] )
  | _ -> ([], [])

let verify g sums =
  let eps = 1.0e-12 in
  let rr, ri = ref_values g.name in
  let rr = Array.of_list (List.map (fun v -> v *. 100.0) rr) in
  let ri = Array.of_list (List.map (fun v -> v *. 100.0) ri) in
  let ok = ref true in
  for i = 1 to g.niter do
    let sr, si = sums.(i) in
    if not (abs_float ((sr -. rr.(i - 1)) /. rr.(i - 1)) < eps) then
      ok := false;
    if not (abs_float ((si -. ri.(i - 1)) /. ri.(i - 1)) < eps) then
      ok := false
  done;
  !ok

(* ---------- main ---------- *)
let () =
  let cls = if Array.length Sys.argv > 1 then Sys.argv.(1) else "S" in
  let g = params_of cls in
  let mask =
    match Sys.getenv_opt "FT_FAR" with
    | None -> 0
    | Some s ->
        match int_of_string_opt s with
        | Some n when n >= 0 && n <= 31 -> n
        | _ -> invalid_arg "FT_FAR must be an integer from 0 to 31"
  in
  let far bit = mask land (1 lsl bit) <> 0 in
  Printf.printf "\n NAS Parallel Benchmarks - FT (OCaml, sequential)\n";
  Printf.printf " Size : %dx%dx%d   Iterations : %d\n" g.nx g.ny g.nz g.niter;
  Printf.printf " FT_FAR : %d (u0=%b u1=%b u2=%b indexmap=%b ex=%b)\n%!"
    mask (far 0) (far 1) (far 2) (far 3) (far 4);
  let t0 = Sys.time () in
  let nmax = max g.nx (max g.ny g.nz) in
  let ur = Array.make nmax 0.0 and ui = Array.make nmax 0.0 in
  let indexmap = compute_indexmap ~far:(far 3) g in
  let ex = compute_exp ~far:(far 4) g in
  let u1 = gen_initial_conditions ~far:(far 1) g in
  let u0 = make_c3 ~far:(far 0) (g.nx * g.ny * g.nz) in
  let u2 = make_c3 ~far:(far 2) (g.nx * g.ny * g.nz) in
  fft_init nmax ur ui;
  fft g ur ui 1 u1 u0;
  (* u0 = FFT(initial conditions) *)
  let sums = Array.make (g.niter + 1) (0.0, 0.0) in
  let loop_start = Unix.gettimeofday () in
  for iter = 1 to g.niter do
    evolve u0 u1 iter indexmap ex;
    (* u1 = u0 * exp(-4 pi^2 alpha t |k|^2) *)
    fft g ur ui (-1) u1 u2;
    (* u2 = IFFT(u1) *)
    checksum g sums iter u2
  done;
  let loop_wall_time = Unix.gettimeofday () -. loop_start in
  let ok = verify g sums in
  let total_cpu_time = Sys.time () -. t0 in
  for iter = 1 to g.niter do
    let sr, si = sums.(iter) in
    Printf.printf "T = %5d     Checksum = %20.12e %20.12e\n" iter sr si
  done;
  Printf.printf
    " Verification: %s\n Total CPU time: %.6f s\n Loop wall time: %.6f s\n%!"
    (if ok then "SUCCESSFUL" else "FAILED")
    total_cpu_time loop_wall_time;
  if not ok then exit 1
