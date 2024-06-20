import { createRoot } from "react-dom/client";
import { createContext, ReactNode, useState } from "react";

// Include the default CSS
import "./Common/base.css";
import SelectDockingFolder from "./Components/SelectDockingFolder";
import { ModelResultsData } from "./Components/DockingResults";

// Create the container
const container = document.createElement("div");
container.id = "root-page";
document.body.appendChild(container);

export type DockingAnalasysContext = {
  selectedFolder: string;
  updateSelectedFolder: (f: string) => void;
  updateCurrentView: (v: ReactNode) => void;
  dockingResults: {
    [modelName: string]: ModelResultsData;
  };
  updateDockingResults: (r: { [modelName: string]: ModelResultsData }) => void;
};

export const DockingContext = createContext<null | DockingAnalasysContext>(
  null
);

function MainView() {
  const [currentView, setCurrentView] = useState<ReactNode>(
    <SelectDockingFolder />
  );
  const [dockingResults, setDockingResults] = useState<{
    [modelName: string]: ModelResultsData;
  } | null>(null);
  const [selectedFolder, setSelectedFolder] = useState("");

  const updateCurrentView = (v: ReactNode) => {
    setCurrentView(v);
  };

  const updateDockingResults = (r: {
    [modelName: string]: ModelResultsData;
  }) => {
    setDockingResults(r);
  };

  const updateSelectedFolder = (f: string) => {
    setSelectedFolder(f);
  };

  return (
    <>
      <DockingContext.Provider
        value={{
          selectedFolder,
          updateSelectedFolder,
          updateCurrentView,
          dockingResults: dockingResults!,
          updateDockingResults,
        }}
      >
        {currentView}
      </DockingContext.Provider>
    </>
  );
}

createRoot(container).render(
  <div className="root-page">
    <MainView />
  </div>
);
