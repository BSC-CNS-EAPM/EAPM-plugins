import { useContext, useEffect, useState } from "react";
import { postData } from "../Common/utils";

import EAPMButton from "./UI/EAPMButton";
import EAPMInput from "./UI/EAPMInput";
import EAPMDiv from "./UI/EAPMDiv";
import { DockingContext } from "..";
import DockingResults from "./DockingResults";

declare global {
  interface Window {
    extensionData: any;
    horus: {
      openExtensionFilePicker: (options: FPOptions) => string;
    };
  }
}

type FPOptions = {
  openFolder?: boolean;
  allowedExtensions?: string[];
  onFileSelect?: (filePath: string) => void;
  onFileConfirm?: (fielPath: string) => void;
};

export default function SelectDockingFolder() {
  const [dockingFolder, setDockingFolder] = useState<string | null>(null);
  const [ligandSeparator, setLigandSeparator] = useState<string>("-");
  const [loading, setLoading] = useState(false);

  const dockingContext = useContext(DockingContext);

  const analyseDocking = async (
    ligandSeparator: string,
    dockingFolder: string
  ) => {
    try {
      setLoading(true);
      const url = window.location.href + "api/get_models";

      const data = await postData(url, { ligandSeparator, dockingFolder });

      if (!data.ok) {
        alert(data.msg);
        return;
      }
      dockingContext?.updateSelectedFolder(dockingFolder);
      dockingContext?.updateDockingResults(data.results);
      dockingContext?.updateCurrentView(<DockingResults />);
    } catch (error) {
      console.error(error);
    } finally {
      setLoading(false);
    }
  };

  // Development
  useEffect(() => {
    if (module.hot) {
      analyseDocking(
        "_",
        "/Users/cdominguez/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Glide_testing/docking"
      );
    }
  }, []);

  useEffect(() => {
    if (window.parent.extensionData) {
      const { blockDockingFolder, blockLigandSeparator } =
        window.parent.extensionData;

      if (blockDockingFolder) {
        setDockingFolder(blockDockingFolder);
        setLigandSeparator(blockLigandSeparator ?? "-");
      }
    }
  }, []);

  const getDockingFolder = () => {
    parent.horus.openExtensionFilePicker({
      openFolder: true,
      onFileConfirm: (filepath: string) => {
        setDockingFolder(filepath);
      },
    });
  };

  return (
    <div className="flex m-auto h-screen text-black">
      <div className="flex flex-col w-fit justify-center m-auto">
        <div className="flex flex-row">
          <EAPMButton className="w-[300px]" onClick={getDockingFolder}>
            Select the Glide results folder
          </EAPMButton>
          <EAPMInput
            placeholder="/path/to/docking/"
            style={{
              width: "500px",
            }}
            value={dockingFolder ?? ""}
            onChange={(e) => {
              setDockingFolder(e.target.value);
            }}
          />
        </div>
        <div className="flex flex-row">
          <EAPMDiv className="w-[300px]">Ligand separator</EAPMDiv>
          <EAPMInput
            placeholder="Ligand separator such as '-' or '_'"
            style={{
              width: "500px",
            }}
            value={ligandSeparator}
            onChange={(e) => {
              setLigandSeparator(e.target.value);
            }}
          />
        </div>
        <EAPMButton
          disabled={loading}
          className={`mt-4 ${
            loading &&
            "bg-gradient-to-r from-orange-500 via-orange-600 to-orange-700 shadow-orange-500/50"
          }`}
          onClick={(_) => {
            analyseDocking(ligandSeparator ?? "-", dockingFolder ?? ".");
          }}
        >
          {loading ? "Analysing..." : "Analyse"}
        </EAPMButton>
      </div>
    </div>
  );
}

export function CurrentSelectedFolder() {
  const dockingContext = useContext(DockingContext);

  return (
    <div className="flex flex-row justify-center gap-2 w-full p-2">
      <EAPMInput
        disabled
        value={dockingContext?.selectedFolder}
        className="w-[80%]"
      />
      <EAPMButton
        onClick={() => {
          dockingContext?.updateCurrentView(<SelectDockingFolder />);
        }}
      >
        Select another folder
      </EAPMButton>
    </div>
  );
}
