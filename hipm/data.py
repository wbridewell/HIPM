# Copyright (c) 2008, Institute for the Study of Learning and Expertise
# All rights reserved.
# For details, see the LICENSE file.

class data:

	def __init__(self, data_file, time_step = None, start_time = 0.0, filter = None, lagramge_flag = False):

		self.file = data_file;
		self.vars = [];
		self.data = {};
		self.time = None;
		self.length = 0;

		if data_file == None: return;
		self.read_from_file(data_file, lagramge_flag = lagramge_flag);

		if time_step <> None:
			self.vars.append("time");
			self.data["time"] = [start_time + i * time_step for i in range(len(self.data[self.vars[0]]))];
			self.time = self.data["time"];
		else:
			self.time = self.data[self.vars[0]];


	def read_from_file(self, data_file, append_flag = False, filter = None, lagramge_flag = False):
		from distutils.text_file import TextFile;
		from string import split, atof;

		f = TextFile(filename = data_file);

		read_vars_flag = not append_flag or (self.length == 0);

		line = f.readline();
		if read_vars_flag:
			if lagramge_flag:
				self.vars = line[-1].split(" ");
			else:
				self.vars = line.split(" ");
			for v in self.vars: self.data[v] = [];

		line_index = -1;
		for line in f.readlines():
			line_index = line_index + 1;
			if (filter <> None) and (line_index not in filter): continue;

			if lagramge_flag:
				vals = line[:-1].split(" ");
			else:
				vals = line.split(" ");
			for i in range(len(vals)): 
				if vals[i] <> "?":
					self.data[self.vars[i]].append(atof(vals[i]));
				else:
					self.data[self.vars[i]].append(vals[i]);

			self.length = self.length + 1;

		f.close();


	def subset(self, filter, file_suffix = ""):

		fdata = data(None);

		fdata.file = self.file[:self.file.rindex(".")] + file_suffix + self.file[self.file.rindex("."):];
		fdata.length = len(filter);

		for v in self.vars:
			fdata.vars.append(v);
			fdata.data[v] = [];
			for i in range(len(self.data[v])):
				if i not in filter: continue;
				fdata.data[v].append(self.data[v][i]);
			if self.data[v] == self.time: fdata.time  = fdata.data[v];
		return fdata;


	def write_to_file(self):

		f = open(self.file, "w");

		first_flag = True;
		for v in self.vars:
			if not first_flag: f.write(" ");
			f.write(v);
			first_flag = False;
		f.write("\n");

		for i in range(self.length):
			first_flag = True;
			for v in self.vars:
				if not first_flag: f.write(" ");
				f.write(str(self.data[v][i]));
				first_flag = False;
			f.write("\n");

		f.close();

	def contains_varname(self, vname): return vname in self.vars;

	def var_index(self, vname):

		for i in range(len(self.vars)):
			if self.vars[i] == vname: return i;
		return -1;

	def value(self, v, t):

		if t >= self.time[-1]: return self.data[v][-1];
		i = 0;
		while self.time[i] <= t: i = i + 1;
		return self.data[v][i - 1];

	def valuei(self, vi, t): return self.value(self.vars[vi], t);


	# degree of fit report for a simulation data set
	def dof_stats(self):
		from math import sqrt;

		def sse(x, y): return sum([(x[i] - y[i]) * (x[i] - y[i]) for i in range(len(x))]);
		def mse(x, y): return sse(x, y) / len(x);
		def rmse(x, y): return sqrt(mse(x, y));

		def remse(x, y):
			def disp(x):
				(sx, sxx) = (0.0, 0.0);
				for i in range(len(x)):
					sx = sx + x[i];
					sxx = sxx + x[i] * x[i];
				return (sxx - sx * sx / len(x)) / (len(x) - 1);
			return mse(x, y) / disp(x);

		def r2(x, y):
			(sx, sy, sxx, sxy, syy) = (0.0, 0.0, 0.0, 0.0, 0.0);
			for i in range(len(x)):
				sx = sx + x[i];
				sy = sy + y[i];
				sxx = sxx + x[i] * x[i];
				syy = syy + y[i] * y[i];
				sxy = sxy + x[i] * y[i];
			ssxx = sxx - sx * sx / len(x);
			ssyy = syy - sy * sy / len(y);
			ssxy = sxy - sx * sy / len(x);
			return (ssxy * ssxy) / (ssxx * ssyy);


		(dall, dall_sim) = ([], []);
		(sum_remse, sum_r2, n) = (0.0, 0.0, 0);
		for v in self.vars[1:]:
			if v[-4:] == "_sim": continue;
			n = n + 1;

			v_sim = v + "_sim";
			[d, d_sim] = [self.data[x] for x in [v, v_sim]];

			dall = dall + d;
			dall_sim = dall_sim + d_sim;

			single_sse = sse(d, d_sim);
			single_mse = mse(d, d_sim);
			single_rmse = rmse(d, d_sim);
			single_remse = remse(d, d_sim);
			single_r2 = r2(d, d_sim);

			sum_remse = sum_remse + single_remse;
			sum_r2 = sum_r2 + single_r2;

			print "%s.SSE = %g" % (v, single_sse);
			print "%s.MSE = %g" % (v, single_mse);
			print "%s.RMSE = %g" % (v, single_rmse);
			print "%s.reMSE = %g" % (v, single_remse);
			print "%s.r2 = %g" % (v, single_r2);
			print "";

		print "AVG.reMSE = %g" % (sum_remse / n);
		print "AVG.r2 = %g" % (sum_r2 / n);
		print "ALL.r2 = %g" % r2(dall, dall_sim);
