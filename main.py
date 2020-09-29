import uproot
import sys

try:
	import ROOT as root
except Exception as ex:
	print("Check if you have ROOT version >= 6.22 installed.")
	raise ex

import numpy as np
import pandas as pd

from config import *


def assemble_DF(tree):
	dataframe = pd.DataFrame({"mJJ":			pd.Series(tree["mJJ"].array()),
							  "deltaYJJ":		pd.Series(tree["deltaYJJ"].array()),
							  "metPt":			pd.Series(tree["metPt"].array()),
							  "ptBalance":		pd.Series(tree["ptBalance"].array()),
							  "subleadJetEta":	pd.Series(tree["subleadJetEta"].array()),
							  "leadJetPt":		pd.Series(tree["leadJetPt"].array()),
							  "photonEta":		pd.Series(tree["photonEta"].array()),
							  "ptBalanceRed":	pd.Series(tree["ptBalanceRed"].array()),
							  "nJets":			pd.Series(tree["nJets"].array()),
							  "sinDeltaPhiJJOver2": pd.Series(tree["sinDeltaPhiJJOver2"].array()),
							  "deltaYJPh":		pd.Series(tree["deltaYJPh"].array()),
							  "weightModified":	pd.Series(tree["weightModified"].array()),
							  "nLeptons":		pd.Series(tree["nLeptons"].array()),
							  "phCentrality":	pd.Series(tree["phCentrality"].array()),})
	
	return dataframe


def apply_filters(dataframe, filter_id=None):
	_len0 = len(dataframe["mJJ"])
	if filter_id == 1:
		dataframe = dataframe[dataframe["nLeptons"] == N_LEP]
		dataframe = dataframe[dataframe["nJets"] > 1]
	elif filter_id == 2:
		dataframe = dataframe[dataframe["nLeptons"] == N_LEP]
		dataframe = dataframe[dataframe["mJJ"] > MJJ_LIM]
		dataframe = dataframe[dataframe["phCentrality"] < PH_CENT]
		dataframe = dataframe[dataframe["nJets"] > 1]
	_len1 = len(dataframe["mJJ"])

	if FILTER == 1:
		filtername = "Z_gamma region"
	elif FILTER == 2:
		filtername = "Signal region"
	else:
		filtername = "_"
	print("Applying {} filter...".format(filtername))
	print("Dropped {} entries after filtering.\n".format(_len0-_len1))
	
	return dataframe.iloc[:, :-2]


def file_read():
	try:
		signal_file = uproot.open(SIGNALPROCFNM)
		bg_file = uproot.open(BGPROCFNM)
	except FileNotFoundError as ex:
		print("Check if you've followed the instrucrions in README.")
		sys.exit()

	signal_trees = [signal_file[TREENM+";1"], signal_file[TREENM+";2"]]
	bg_trees = [bg_file[TREENM+";10"], bg_file[TREENM+";11"]]

	raw_bg_data = pd.concat([assemble_DF(bg_trees[0]), assemble_DF(bg_trees[1])], ignore_index=True)
	raw_signal_data = pd.concat([assemble_DF(signal_trees[0]), assemble_DF(signal_trees[1])], ignore_index=True)

	filtered_bg_data = apply_filters(raw_bg_data, filter_id=FILTER)
	filtered_signal_data = apply_filters(raw_signal_data, filter_id=FILTER)

	return filtered_bg_data, filtered_signal_data



def setup_layout(canvas, bg_hist, signal_hist):
	canvas.SetTicks(1, 1)					# добавляет метки на верхнюю и правую оси
	
	bg_hist.SetStats(False)					# убирает сводную таблицу
	signal_hist.SetStats(False)
	bg_hist.SetLineWidth(2)	
	bg_hist.SetLineColor(2)
	bg_hist.SetFillColor(2)
	bg_hist.SetFillStyle(3004)				# устанавливает стиль заливки

	signal_hist.SetLineWidth(2)
	signal_hist.SetLineColor(4)
	signal_hist.SetFillColorAlpha(4, 0.35)	# добавляет прозрачность заливки

	legend=root.TLegend(0.72, 0.77, 0.9, 0.9)	# создаёт объект легенды
	legend.AddEntry(signal_hist,"Signal","f")
	legend.AddEntry(bg_hist,"Background","f")

	latex = root.TLatex()
	latex.SetNDC()
	latex.SetTextSize(0.035)

	return canvas, bg_hist, signal_hist, legend, latex


def plot(filtered_bg_data, filtered_signal_data):
	for key in filtered_signal_data.keys():
		bg_column = filtered_bg_data[key]
		signal_column = filtered_signal_data[key]
		print("BG: {} entries.".format(len(bg_column)))
		print("SIGN: {} entries.".format(len(signal_column)))
		print("Total: {} entries.".format(len(bg_column) + len(signal_column)))

		xlim = np.max([np.max(bg_column), np.max(signal_column)])
		xzero = 0
		bins = DEF_BINS
		xlabel = ""

		if key in ("subleadJetEta", "photonEta"):
			xzero = -xlim
		elif key == "nJets":
			bins = 8

		if key in ("mJJ", "metPt", "leadJetPt"):
			xlabel = "GeV"

		canvas = root.TCanvas("canvas", "CANVAS")
		canvas.SetTitle(key)
		bg_hist = root.TH1F("", key, bins, xzero, xlim)
		bg_hist.GetXaxis().SetTitle(xlabel)
		signal_hist = root.TH1F("", key, bins, xzero, xlim)
		
		canvas, bg_hist, signal_hist, legend, latex = setup_layout(canvas, bg_hist, signal_hist)


		for bg_value, weight in zip(bg_column, filtered_bg_data["weightModified"]):
			bg_hist.Fill(bg_value)
		for signl_value, weight in zip(signal_column, filtered_signal_data["weightModified"]):
			signal_hist.Fill(signl_value)

		if FILTER == 1:
			text = "Z_{#gamma} region"

			if key in ("mJJ", "deltaYJJ", "metPt", "subleadJetEta", "leadJetPt", "sinDeltaPhiJJOver2", "deltaYJPh"):
				bg_hist.DrawNormalized()
				signal_hist.DrawNormalized("same")
			else:
				signal_hist.DrawNormalized()
				bg_hist.DrawNormalized("same")	
		elif FILTER == 2:
			text = "Signal region"
			if key in ("mJJ", "deltaYJJ", "metPt", "photonEta", "sinDeltaPhiJJOver2", "deltaYJPh"):
				bg_hist.DrawNormalized()
				signal_hist.DrawNormalized("same")
			else:
				signal_hist.DrawNormalized()
				bg_hist.DrawNormalized("same")	
		else:
			text = ""

		# bg_hist.DrawNormalized()
		# signal_hist.DrawNormalized("same")

		legend.Draw()
		latex.DrawLatex(0.75, 0.73, text)
		canvas.Update()
		t = input()	

		bg_hist = None
		signal_hist = None
		canvas = None


if __name__ == "__main__":
	plot(*file_read())
