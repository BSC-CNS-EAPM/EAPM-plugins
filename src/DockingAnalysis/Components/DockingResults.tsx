import { useContext, useState } from "react";
import { DockingContext } from "..";
import EAPMButton from "./UI/EAPMButton";
import { CurrentSelectedFolder } from "./SelectDockingFolder";
import { EAPMTable } from "./UI/EAPMTable";

import Modal from "react-modal";
import { LazyLog } from "@melloware/react-logviewer";
import EAPMDiv from "./UI/EAPMDiv";

export type ModelResultsData = {
  [ligandName: string]: LigandResultsData;
};

export type LigandResultsData = {
  ligandName: string;
  separator: string;
  resultsCSVPath: string;
  logPath: string;
  logSubjobsPath: string;
  csvData: DockingResultsCSV;
};

type DockingResultsCSV = {
  columns: string[];
  rows: any[];
};

export default function DockingResults() {
  const [currentModel, setCurrentModel] = useState<string | null>(null);
  const [currentLigand, setCurrentLigand] = useState<string | null>(null);
  const dockingContext = useContext(DockingContext);
  const results = dockingContext?.dockingResults || {};
  const models = Object.keys(results);
  const ligands = currentModel ? Object.keys(results[currentModel] || {}) : [];

  const handleModelClick = (model: string) => {
    setCurrentModel(model);
    setCurrentLigand(null); // Reset current ligand when changing model
  };

  const handleLigandClick = (ligand: string) => {
    setCurrentLigand(ligand);
  };

  return (
    <>
      <CurrentSelectedFolder />
      <div className="flex flex-col p-4 items-start justify-start w-full">
        <div className="flex flex-row gap-1 items-center">
          <EAPMDiv>Models</EAPMDiv>
          {models.map((model) => (
            <EAPMButton
              key={model}
              setClass={false}
              onClick={() => handleModelClick(model)}
              className={`${
                currentModel === model
                  ? "bg-blue-500 text-white"
                  : "border-blue-500 text-blue-500"
              }`}
            >
              {model}
            </EAPMButton>
          ))}
        </div>
        {currentModel && (
          <div className="flex flex-row gap-1 items-center">
            <EAPMDiv>Ligands</EAPMDiv>
            {ligands.map((ligand) => (
              <EAPMButton
                setClass={false}
                key={ligand}
                onClick={() => handleLigandClick(ligand)}
                className={`${
                  currentLigand === ligand
                    ? "bg-blue-500 text-white"
                    : "border-blue-500 text-blue-500"
                }`}
              >
                {ligand}
              </EAPMButton>
            ))}
          </div>
        )}
        {currentModel && currentLigand && (
          <SelectedLigandView
            model={currentModel}
            ligand={results[currentModel]![currentLigand]!}
          />
        )}
      </div>
    </>
  );
}
type SelectedLigandViewProps = {
  model: string;
  ligand: LigandResultsData;
};

function SelectedLigandView(props: SelectedLigandViewProps) {
  const [logsOpen, setLogsOpen] = useState(false);

  const logsURL = () => {
    return (
      window.location.href +
      "api/get_file_contents?path=" +
      encodeURIComponent(props.ligand.logPath)
    );
  };

  return (
    <>
      <div className="flex flex-row gap-1">
        <EAPMButton
          buttonColor="orange"
          className="flex flex-row gap-1 items-center"
          onClick={(_) => {
            setLogsOpen(true);
          }}
        >
          <LogFileIcon />
          View logs
        </EAPMButton>
        <EAPMButton
          buttonColor="orange"
          className="flex flex-row gap-1 items-center"
          onClick={(_) => {
            setLogsOpen(true);
          }}
        >
          <ExtractIcon />
          Extract all poses
        </EAPMButton>
        <EAPMButton
          buttonColor="orange"
          className="flex flex-row gap-1 items-center"
          onClick={(_) => {
            setLogsOpen(true);
          }}
        >
          <AnalyseIcon />
          Analyse
        </EAPMButton>
      </div>
      <EAPMTable results={props.ligand.csvData} />
      <Modal
        style={{
          position: "absolute",
          zIndex: Number.MAX_SAFE_INTEGER,
        }}
        isOpen={logsOpen}
        onRequestClose={() => {
          setLogsOpen(false);
        }}
        contentLabel="Ligand logs"
      >
        <LazyLog url={logsURL()} enableSearchNavigation enableSearch />
      </Modal>
    </>
  );
}

function LogFileIcon() {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      fill="none"
      viewBox="0 0 24 24"
      strokeWidth={1.5}
      stroke="currentColor"
      className="size-6"
    >
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        d="M19.5 14.25v-2.625a3.375 3.375 0 0 0-3.375-3.375h-1.5A1.125 1.125 0 0 1 13.5 7.125v-1.5a3.375 3.375 0 0 0-3.375-3.375H8.25m0 12.75h7.5m-7.5 3H12M10.5 2.25H5.625c-.621 0-1.125.504-1.125 1.125v17.25c0 .621.504 1.125 1.125 1.125h12.75c.621 0 1.125-.504 1.125-1.125V11.25a9 9 0 0 0-9-9Z"
      />
    </svg>
  );
}

function ExtractIcon() {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      fill="none"
      viewBox="0 0 24 24"
      strokeWidth={1.5}
      stroke="currentColor"
      className="size-6"
    >
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        d="M13.5 6H5.25A2.25 2.25 0 0 0 3 8.25v10.5A2.25 2.25 0 0 0 5.25 21h10.5A2.25 2.25 0 0 0 18 18.75V10.5m-10.5 6L21 3m0 0h-5.25M21 3v5.25"
      />
    </svg>
  );
}

function AnalyseIcon() {
  return (
    <svg
      xmlns="http://www.w3.org/2000/svg"
      fill="none"
      viewBox="0 0 24 24"
      strokeWidth={1.5}
      stroke="currentColor"
      className="size-6"
    >
      <path
        strokeLinecap="round"
        strokeLinejoin="round"
        d="m21 21-5.197-5.197m0 0A7.5 7.5 0 1 0 5.196 5.196a7.5 7.5 0 0 0 10.607 10.607Z"
      />
    </svg>
  );
}
